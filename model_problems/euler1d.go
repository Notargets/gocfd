package model_problems

import (
	"fmt"
	"math"
	"sync"
	"time"

	"github.com/notargets/avs/chart2d"
	utils2 "github.com/notargets/avs/utils"

	"github.com/notargets/gocfd/DG1D"
	"github.com/notargets/gocfd/utils"
)

type Euler1D struct {
	// Input parameters
	CFL, FinalTime  float64
	El              *DG1D.Elements1D
	RHSOnce         sync.Once
	State           *FieldState
	Rho, RhoU, Ener utils.Matrix
	In, Out         *State
	plotOnce        sync.Once
	chart           *chart2d.Chart2D
	colorMap        *utils2.ColorMap
	frameCount      int
}

type FieldState struct {
	Gamma                    float64
	U, Q, Pres, CVel, LM, Ht utils.Matrix
	Temp                     utils.Matrix
}

func NewFieldState() (fs *FieldState) {
	return &FieldState{}
}

func (fs *FieldState) Update(Rho, RhoU, Ener utils.Matrix) {
	/*
		Calorically Perfect Gas, R = 1
	*/
	var (
		Gamma = 1.4
		Cv    = 1. / (Gamma - 1.)
		Cp    = Gamma * Cv
	)
	fs.U = RhoU.Copy().ElDiv(Rho)                   // Velocity
	fs.Q = fs.U.Copy().POW(2).ElMul(Rho).Scale(0.5) // 1/2 * Rho * U^2
	// (gamma-1.0)*(Ener - 0.5*(rhou).^2./rho)
	fs.Pres = Ener.Copy().Subtract(fs.Q).Scale(fs.Gamma - 1.) // Pressure
	// sqrt(gamma*pres./rho)
	fs.CVel = fs.Pres.Copy().ElDiv(Rho).Scale(fs.Gamma).Apply(math.Sqrt)
	//  abs(rhou./rho)+cvel
	fs.LM = fs.U.Copy().Apply(math.Abs).Add(fs.CVel)
	//	Temp = (Ener - 0.5*(rhou).^2./rho)./rho
	fs.Temp = Ener.Copy().Subtract(fs.Q).ElDiv(Rho)
	fs.Ht = fs.Q.Copy().ElDiv(Rho).Add(fs.Temp.Copy().Scale(Cp))
	fs.U.SetReadOnly("U")
	fs.Q.SetReadOnly("Q")
	fs.Pres.SetReadOnly("Pres")
	fs.CVel.SetReadOnly("CVel")
	fs.LM.SetReadOnly("LM")
	fs.Ht.SetReadOnly("Ht")
}

func (fs *FieldState) Print() {
	fmt.Println(fs.U.Print("U"))
	fmt.Println(fs.Q.Print("Q"))
	fmt.Println(fs.Pres.Print("Pres"))
	fmt.Println(fs.CVel.Print("CVel"))
	fmt.Println(fs.LM.Print("LM"))
	fmt.Println(fs.Ht.Print("Ht"))
}

func NewEuler1D(CFL, FinalTime, XMax float64, N, K int) (c *Euler1D) {
	VX, EToV := DG1D.SimpleMesh1D(0, XMax, K)
	c = &Euler1D{
		CFL:       CFL,
		State:     NewFieldState(),
		FinalTime: FinalTime,
		El:        DG1D.NewElements1D(N, VX, EToV),
	}
	c.State.Gamma = 1.4
	c.In = NewStateP(c.State.Gamma, 1, 0, 1)
	c.Out = NewStateP(c.State.Gamma, 0.125, 0, 0.1)
	fmt.Printf("Euler Equations in 1 Dimension\nSolving Sod's Shock Tube\n")
	fmt.Printf("CFL = %8.4f, Polynomial Degree N = %d (1 is linear), Num Elements K = %d\n\n\n", CFL, N, K)
	prob := "SOD"
	switch prob {
	case "SOD":
		c.InitializeSOD()
	case "FS":
		c.Out = c.In
		c.InitializeFS()
	case "Collision":
		fallthrough
	default:
		c.Out = c.In
		c.InitializeSOD()
	}
	return
}

func (c *Euler1D) InitializeFS() {
	var (
		el = c.El
		FS = c.In
	)
	c.Rho = utils.NewMatrix(el.Np, el.K).AddScalar(FS.Rho)
	c.RhoU = utils.NewMatrix(el.Np, el.K).AddScalar(FS.RhoU)
	c.Ener = utils.NewMatrix(el.Np, el.K).AddScalar(FS.Ener)
}

func (c *Euler1D) InitializeSOD() {
	var (
		/*
			In                = NewStateP(c.Gamma, 1, 0, 1)
			Out               = NewStateP(c.Gamma, 0.125, 0, 0.1)
		*/
		el                = c.El
		MassMatrix, VtInv utils.Matrix
		err               error
		npOnes            = utils.NewVectorConstant(el.Np, 1)
		s                 = c.State
	)
	if VtInv, err = el.V.Transpose().Inverse(); err != nil {
		panic(err)
	}
	MassMatrix = VtInv.Mul(el.Vinv)
	CellCenterRValues := MassMatrix.Mul(el.X).SumCols().Scale(0.5)
	cx := npOnes.Outer(CellCenterRValues)
	leftHalf := cx.Find(utils.Less, 0.5, false)
	rightHalf := cx.Find(utils.GreaterOrEqual, 0.5, false)
	// Initialize field variables
	c.Rho = utils.NewMatrix(el.Np, el.K)
	c.Rho.AssignScalar(leftHalf, 1)
	c.Rho.AssignScalar(rightHalf, 0.125)
	c.RhoU = utils.NewMatrix(el.Np, el.K)
	c.RhoU.Scale(0)
	rDiv := 1. / (s.Gamma - 1.)
	c.Ener = utils.NewMatrix(el.Np, el.K)
	c.Ener.AssignScalar(leftHalf, rDiv)
	c.Ener.AssignScalar(rightHalf, 0.1*rDiv)
}

func (c *Euler1D) Run(showGraph bool, graphDelay ...time.Duration) {
	var (
		el           = c.El
		logFrequency = 50
		//s             = c.State
	)
	xmin := el.X.Row(1).Subtract(el.X.Row(0)).Apply(math.Abs).Min()
	var Time, dt float64
	var tstep int
	for Time < c.FinalTime {
		/*
			Third Order Runge-Kutta time advancement
		*/
		// SSP RK Stage 1
		rhsRho, rhsRhoU, rhsEner := c.RHS(&c.Rho, &c.RhoU, &c.Ener)
		c.Plot(Time, showGraph, graphDelay)
		dt = c.CalculateDT(xmin, Time)
		rho1 := c.Rho.Copy().Add(rhsRho.Scale(dt))
		rhou1 := c.RhoU.Copy().Add(rhsRhoU.Scale(dt))
		ener1 := c.Ener.Copy().Add(rhsEner.Scale(dt))

		// SSP RK Stage 2
		rhsRho, rhsRhoU, rhsEner = c.RHS(&rho1, &rhou1, &ener1)
		rho2 := c.Rho.Copy().Scale(3.).Add(rho1).Add(rhsRho.Scale(dt)).Scale(0.25)
		rhou2 := c.RhoU.Copy().Scale(3.).Add(rhou1).Add(rhsRhoU.Scale(dt)).Scale(0.25)
		ener2 := c.Ener.Copy().Scale(3.).Add(ener1).Add(rhsEner.Scale(dt)).Scale(0.25)

		// SSP RK Stage 3
		rhsRho, rhsRhoU, rhsEner = c.RHS(&rho2, &rhou2, &ener2)
		c.Rho.Add(rho2.Scale(2)).Add(rhsRho.Scale(2. * dt)).Scale(1. / 3.)
		c.RhoU.Add(rhou2.Scale(2)).Add(rhsRhoU.Scale(2. * dt)).Scale(1. / 3.)
		c.Ener.Add(ener2.Scale(2)).Add(rhsEner.Copy().Scale(2. * dt)).Scale(1. / 3.)

		Time += dt
		tstep++
		isDone := math.Abs(Time-c.FinalTime) < 0.00001
		if tstep%logFrequency == 0 || isDone {
			fmt.Printf("Time = %8.4f, max_resid[%d] = %8.4f, emin = %8.6f, emax = %8.6f\n", Time, tstep, rhsEner.Max(), c.Ener.Min(), c.Ener.Max())
			if isDone {
				for {
					time.Sleep(time.Second)
				}
			}
		}
	}
	return
}

func (c *Euler1D) CalculateDT(xmin, Time float64) (dt float64) {
	var (
		s = c.State
	)
	// min(xmin ./ (abs(U) +C))
	Factor := s.U.Copy().Apply(math.Abs).Add(s.CVel).Apply(func(val float64) float64 { return xmin / val })
	dt = c.CFL * Factor.Min()
	if dt+Time > c.FinalTime {
		dt = c.FinalTime - Time
	}
	return
}

func (c *Euler1D) RHS(Rhop, RhoUp, Enerp *utils.Matrix) (rhsRho, rhsRhoU, rhsEner utils.Matrix) {
	var (
		el                                                 = c.El
		nrF, ncF                                           = el.Nfp * el.NFaces, el.K
		s                                                  = c.State
		dRho, dRhoU, dEner, dRhoF, dRhoUF, dEnerF, LFcDiv2 utils.Matrix
		Rho, RhoU, Ener                                    = *Rhop, *RhoUp, *Enerp
		RhoF, RhoUF, EnerF                                 utils.Matrix
		limiter                                            = true
		slopeLimiterM                                      = 20.
	)
	if limiter {
		// Slope Limit the solution fields
		*Rhop = el.SlopeLimitN(*Rhop, slopeLimiterM)
		*RhoUp = el.SlopeLimitN(*RhoUp, slopeLimiterM)
		*Enerp = el.SlopeLimitN(*Enerp, slopeLimiterM)
		Rho, RhoU, Ener = *Rhop, *RhoUp, *Enerp
	}
	s.Update(Rho, RhoU, Ener)
	RhoF = RhoU.Copy()
	RhoUF = s.Q.Copy().Scale(2.).Add(s.Pres)
	EnerF = Ener.Copy().Add(s.Pres).ElMul(s.U)
	// Face jumps in primary and flux variables
	fJump := func(U utils.Matrix) (dU utils.Matrix) {
		dU = U.Subset(el.VmapM, nrF, ncF).Subtract(U.Subset(el.VmapP, nrF, ncF))
		return
	}
	dRho = fJump(Rho)
	dRhoU = fJump(RhoU)
	dEner = fJump(Ener)
	dRhoF = fJump(RhoF)
	dRhoUF = fJump(RhoUF)
	dEnerF = fJump(EnerF)
	// Lax-Friedrichs flux component is always used divided by 2, so we pre-scale it
	LFcDiv2 = s.LM.Subset(el.VmapM, nrF, ncF).Apply2(s.LM.Subset(el.VmapP, nrF, ncF), math.Max).Scale(0.5)

	// Compute fluxes at interfaces
	dRhoF.ElMul(el.NX).Scale(0.5).Subtract(LFcDiv2.Copy().ElMul(dRho))
	dRhoUF.ElMul(el.NX).Scale(0.5).Subtract(LFcDiv2.Copy().ElMul(dRhoU))
	dEnerF.ElMul(el.NX).Scale(0.5).Subtract(LFcDiv2.ElMul(dEner))

	c.BoundaryConditions(Rho, RhoU, Ener, RhoF, RhoUF, EnerF, &dRhoF, &dRhoUF, &dEnerF)

	// RHS Computation
	rhsRho = el.Rx.Copy().Scale(-1.).ElMul(el.Dr.Mul(RhoF)).Add(el.LIFT.Mul(dRhoF.ElMul(el.FScale)))
	rhsRhoU = el.Rx.Copy().Scale(-1.).ElMul(el.Dr.Mul(RhoUF)).Add(el.LIFT.Mul(dRhoUF.ElMul(el.FScale)))
	rhsEner = el.Rx.Copy().Scale(-1.).ElMul(el.Dr.Mul(EnerF)).Add(el.LIFT.Mul(dEnerF.ElMul(el.FScale)))
	return
}

func (c *Euler1D) BoundaryConditions(Rho, RhoU, Ener, RhoF, RhoUF, EnerF utils.Matrix, dRhoF, dRhoUF, dEnerF *utils.Matrix) {
	var (
		s  = c.State
		el = c.El
		// Sod's problem: Shock tube with jump in middle
		In  = c.In
		Out = c.Out
	)

	// Boundary conditions for Sod's problem
	// Inflow
	lmI := s.LM.SubsetVector(el.VmapI).Scale(0.5)
	nxI := el.NX.SubsetVector(el.MapI)
	bFunc(dRhoF, Rho, RhoF, lmI, nxI, In.Rho, In.RhoF, el.MapI, el.VmapI)
	bFunc(dRhoUF, RhoU, RhoUF, lmI, nxI, In.RhoU, In.RhoUF, el.MapI, el.VmapI)
	bFunc(dEnerF, Ener, EnerF, lmI, nxI, In.Ener, In.EnerF, el.MapI, el.VmapI)

	// Outflow
	lmO := s.LM.SubsetVector(el.VmapO).Scale(0.5)
	nxO := el.NX.SubsetVector(el.MapO)
	bFunc(dRhoF, Rho, RhoF, lmO, nxO, Out.Rho, Out.RhoF, el.MapO, el.VmapO)
	bFunc(dRhoUF, RhoU, RhoUF, lmO, nxO, Out.RhoU, Out.RhoUF, el.MapO, el.VmapO)
	bFunc(dEnerF, Ener, EnerF, lmO, nxO, Out.Ener, Out.EnerF, el.MapO, el.VmapO)
}

func bFunc(dUF *utils.Matrix, U, UF utils.Matrix, lm, nx utils.Vector, uIO, ufIO float64, mapi, vmap utils.Index) {
	// Characteristic BC using freestream conditions at boundary
	var (
		dUFVec utils.Matrix
	)
	dUFVec = nx.Outer(UF.SubsetVector(vmap).Subtract(utils.NewVectorConstant(len(vmap), ufIO))).Scale(0.5)
	dUFVec.Subtract(lm.Outer(U.SubsetVector(vmap).Subtract(utils.NewVectorConstant(len(vmap), uIO))))
	dUF.AssignVector(mapi, dUFVec)
	return
}

type State struct {
	Gamma, Rho, RhoU, Ener float64
	RhoF, RhoUF, EnerF     float64
}

func NewState(gamma, rho, rhoU, ener float64) (s *State) {
	p := ener * (gamma - 1.)
	q := 0.5 * utils.POW(rhoU, 2) / rho
	u := rhoU / rho
	return &State{
		Gamma: gamma,
		Rho:   rho,
		RhoU:  rhoU,
		Ener:  ener,
		RhoF:  rhoU,
		RhoUF: 2*q + p,
		EnerF: (ener + q + p) * u,
	}
}

func NewStateP(gamma, rho, rhoU, p float64) *State {
	ener := p / (gamma - 1.)
	return NewState(gamma, rho, rhoU, ener)
}

func (s *State) Print() (o string) {
	o = fmt.Sprintf("Rho = %v\nP = %v\nE = %v\nRhoU = %v\nRhoUF = %v\n",
		s.Rho, s.Ener*(s.Gamma-1.), s.Ener, s.RhoU, s.RhoUF)
	return
}

func (c *Euler1D) Plot(timeT float64, showGraph bool, graphDelay []time.Duration) {
	var (
		el         = c.El
		fmin, fmax = float32(-0.1), float32(2.6)
	)
	if !showGraph {
		return
	}
	c.plotOnce.Do(func() {
		c.chart = chart2d.NewChart2D(1920, 1280, float32(el.X.Min()), float32(el.X.Max()), fmin, fmax)
		c.colorMap = utils2.NewColorMap(-1, 1, 1)
		go c.chart.Plot()
	})
	pSeries := func(field utils.Matrix, name string, color float32) {
		if err := c.chart.AddSeries(name, el.X.Transpose().RawMatrix().Data, field.Transpose().RawMatrix().Data,
			chart2d.NoGlyph, chart2d.Solid, c.colorMap.GetRGB(color)); err != nil {
			panic("unable to add graph series")
		}
	}
	pSeries(c.Rho, "Rho", -0.7)
	pSeries(c.RhoU, "RhoU", 0.0)
	pSeries(c.Ener, "Ener", 0.7)
	//pSeries(c.State.U, "U", 0.8)
	//pSeries(c.State.Temp, "Temp", 0.9)
	c.frameCount++
	check := int(math.Log10(float64(el.K * el.Np / 5)))
	if c.frameCount%check == 0 {
		AddAnalyticSod(c.chart, c.colorMap, timeT)
	}
	if len(graphDelay) != 0 {
		time.Sleep(graphDelay[0])
	}
}

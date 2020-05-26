package model_problems

import (
	"fmt"
	"math"
	"sync"
	"time"

	"github.com/notargets/gocfd/sod_shock_tube"

	"github.com/notargets/avs/chart2d"
	utils2 "github.com/notargets/avs/utils"

	"github.com/notargets/gocfd/DG1D"
	"github.com/notargets/gocfd/utils"
)

type EulerDFR struct {
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
	fluxType        string
	frameCount      int
}

func NewEulerDFR(CFL, FinalTime, XMax float64, N, K int) (c *EulerDFR) {
	VX, EToV := DG1D.SimpleMesh1D(0, XMax, K)
	c = &EulerDFR{
		CFL:       CFL,
		State:     NewFieldState(),
		FinalTime: FinalTime,
		El:        DG1D.NewElements1D(N, VX, EToV),
	}
	c.State.Gamma = 1.4
	c.In = NewStateP(c.State.Gamma, 1, 0, 1)
	c.Out = NewStateP(c.State.Gamma, 0.125, 0, 0.1)
	c.fluxType = "LF"
	c.fluxType = "Roe"
	fmt.Printf("Euler Equations in 1 Dimension\nSolving Sod's Shock Tube\nDirect Flux Reconstruction, %s Flux\n", c.fluxType)
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

func (c *EulerDFR) InitializeFS() {
	var (
		el = c.El
		FS = c.In
	)
	c.Rho = utils.NewMatrix(el.Np, el.K).AddScalar(FS.Rho)
	c.RhoU = utils.NewMatrix(el.Np, el.K).AddScalar(FS.RhoU)
	c.Ener = utils.NewMatrix(el.Np, el.K).AddScalar(FS.Ener)
}

func (c *EulerDFR) InitializeSOD() {
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

func (c *EulerDFR) Run(showGraph bool, graphDelay ...time.Duration) {
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

func (c *EulerDFR) CalculateDT(xmin, Time float64) (dt float64) {
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

func (c *EulerDFR) RHS(Rhop, RhoUp, Enerp *utils.Matrix) (rhsRho, rhsRhoU, rhsEner utils.Matrix) {
	var (
		el                 = c.El
		s                  = c.State
		fRho, fRhoU, fEner utils.Matrix
		Rho, RhoU, Ener    = *Rhop, *RhoUp, *Enerp
		RhoF, RhoUF, EnerF utils.Matrix
		limiter            = true
		slopeLimiterM      = 20.
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

	switch c.fluxType {
	case "LF":
		fRho, fRhoU, fEner = c.LFFlux(Rho, RhoU, Ener, RhoF, RhoUF, EnerF)
	case "Roe":
		fRho, fRhoU, fEner = c.RoeFlux(Rho, RhoU, Ener, RhoF, RhoUF, EnerF)
	}

	// Set face flux within global flux
	RhoF.AssignVector(el.VmapM, fRho)
	RhoUF.AssignVector(el.VmapM, fRhoU)
	EnerF.AssignVector(el.VmapM, fEner)

	c.BoundaryConditions(Rho, RhoU, Ener, RhoF, RhoUF, EnerF, &fRho, &fRhoU, &fEner)

	// Calculate RHS
	rhsRho = el.Dr.Mul(RhoF).Scale(-1).ElMul(el.Rx)
	rhsRhoU = el.Dr.Mul(RhoUF).Scale(-1).ElMul(el.Rx)
	rhsEner = el.Dr.Mul(EnerF).Scale(-1).ElMul(el.Rx)
	return
}

func (c *EulerDFR) BoundaryConditions(Rho, RhoU, Ener, RhoF, RhoUF, EnerF utils.Matrix, dRhoF, dRhoUF, dEnerF *utils.Matrix) {
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

func (c *EulerDFR) Plot(timeT float64, showGraph bool, graphDelay []time.Duration) {
	var (
		el         = c.El
		fmin, fmax = float32(-1.5), float32(5)
	)
	if !showGraph {
		return
	}
	c.plotOnce.Do(func() {
		c.chart = chart2d.NewChart2D(1920, 1280, float32(el.X.Min()), float32(el.X.Max()), fmin, fmax)
		c.colorMap = utils2.NewColorMap(-1, 1, 1)
		go c.chart.Plot()
	})
	pSeries := func(field utils.Matrix, name string, color float32, gl chart2d.GlyphType) {
		if err := c.chart.AddSeries(name, el.X.Transpose().RawMatrix().Data, field.Transpose().RawMatrix().Data,
			gl, chart2d.Solid, c.colorMap.GetRGB(color)); err != nil {
			panic("unable to add graph series")
		}
	}
	pSeries(c.Rho, "Rho", -0.7, chart2d.NoGlyph)
	pSeries(c.RhoU, "RhoU", 0.0, chart2d.NoGlyph)
	pSeries(c.Ener, "Ener", 0.6, chart2d.NoGlyph)
	pSeries(c.State.U, "U", 0.8, chart2d.NoGlyph)
	pSeries(c.State.Temp, "Temp", 0.9, chart2d.NoGlyph)
	//pSeries(c.State.Ht, "Ht", -0.9, chart2d.CrossGlyph)
	c.frameCount++
	check := int(math.Log10(float64(el.K * el.Np / 5)))
	if c.frameCount%check == 0 {
		AddAnalyticSod(c.chart, c.colorMap, timeT)
	}
	if len(graphDelay) != 0 {
		time.Sleep(graphDelay[0])
	}
}

func AddAnalyticSod(chart *chart2d.Chart2D, colorMap *utils2.ColorMap, timeT float64) {
	X, Rho, P, U, E := sod_shock_tube.SOD_calc(timeT)
	_, _ = P, U
	if err := chart.AddSeries("ExactRho", X, Rho, chart2d.XGlyph, chart2d.NoLine, colorMap.GetRGB(0.8)); err != nil {
		panic("unable to add exact solution Rho")
	}
	if err := chart.AddSeries("ExactE", X, E, chart2d.XGlyph, chart2d.NoLine, colorMap.GetRGB(1.0)); err != nil {
		panic("unable to add exact solution E")
	}
}

func (c *EulerDFR) LFFlux(Rho, RhoU, Ener, RhoF, RhoUF, EnerF utils.Matrix) (fRho, fRhoU, fEner utils.Matrix) {
	var (
		el       = c.El
		s        = c.State
		nrF, ncF = el.Nfp * el.NFaces, el.K
	)
	// Compute Lax-Friedrichs flux
	// Face jumps in primary and flux variables
	fJump := func(U utils.Matrix) (dU utils.Matrix) {
		dU = U.Subset(el.VmapM, nrF, ncF).Subtract(U.Subset(el.VmapP, nrF, ncF)).ElMul(el.NX)
		return
	}
	// Face flux average
	fAve := func(U utils.Matrix) (Uavg utils.Matrix) {
		Uavg = U.Subset(el.VmapM, nrF, ncF).Add(U.Subset(el.VmapP, nrF, ncF)).Scale(0.5)
		return
	}
	// Max eigenvalue
	LFc := s.LM.Subset(el.VmapM, nrF, ncF).Apply2(s.LM.Subset(el.VmapP, nrF, ncF), math.Max)
	// Compute numerical flux at faces
	fRho = fAve(RhoF).Add(fJump(Rho).ElMul(LFc).Scale(0.5))
	fRhoU = fAve(RhoUF).Add(fJump(RhoU).ElMul(LFc).Scale(0.5))
	fEner = fAve(EnerF).Add(fJump(Ener).ElMul(LFc).Scale(0.5))
	return
}
func (c *EulerDFR) RoeFlux(Rho, RhoU, Ener, RhoF, RhoUF, EnerF utils.Matrix) (fRho, fRhoU, fEner utils.Matrix) {
	var (
		el       = c.El
		nrF, ncF = el.Nfp * el.NFaces, el.K
		s        = c.State
	)
	fL := func(U utils.Matrix) (Ul utils.Matrix) {
		Ul = U.Subset(el.VmapM, nrF, ncF)
		return
	}
	fR := func(U utils.Matrix) (Ul utils.Matrix) {
		Ul = U.Subset(el.VmapP, nrF, ncF)
		return
	}
	// Face jumps in primary and flux variables
	fJump := func(U utils.Matrix) (dU utils.Matrix) {
		dU = U.Subset(el.VmapP, nrF, ncF).Subtract(U.Subset(el.VmapM, nrF, ncF))
		return
	}
	// Face average
	fAve := func(U utils.Matrix) (Uavg utils.Matrix) {
		Uavg = U.Subset(el.VmapP, nrF, ncF).Add(U.Subset(el.VmapM, nrF, ncF)).Scale(0.5)
		return
	}
	/*
		Calculate the Roe Averaged variables
	*/
	RhoL, RhoR := fL(Rho), fR(Rho)
	UL, UR := fL(s.U), fR(s.U)
	HtL, HtR := fL(s.Ht), fR(s.Ht)
	RhoRL := RhoL.Copy().Apply2(RhoR, func(rhol, rhor float64) (res float64) {
		res = math.Sqrt(rhol * rhor)
		return
	})
	roeAve := func(uL, uR utils.Matrix) (uRL utils.Matrix) {
		uRL = RhoL.Copy().Apply4(uL, RhoR, uR, func(rl, ul, rr, ur float64) (res float64) {
			var (
				srl, srr = math.Sqrt(rl), math.Sqrt(rr)
			)
			res = (srl*ul + srr*ur) / (srl + srr)
			return
		})
		return
	}
	URL := roeAve(UL, UR)
	HtRL := roeAve(HtL, HtR)
	aRL := HtRL.Copy().Apply2(URL, func(htrl, url float64) (res float64) {
		res = math.Sqrt((s.Gamma - 1) * (htrl - url*url*0.5))
		return
	})
	DelRho := fJump(Rho)
	DelU := fJump(s.U)
	DelP := fJump(s.Pres)
	// Phi is the Harten entropy correction - it modifies the eigenvalues to eliminate aphysical solutions
	phi := func(eig, del float64) (res float64) {
		absLam := math.Abs(eig)
		if absLam > del {
			res = absLam
		} else {
			res = (eig*eig + del*del) / (2 * del)
		}
		return
	}
	eRho := URL.Copy().Apply8(aRL, RhoRL, HtRL, DelU, DelRho, DelP, el.NX, func(url, arl, rhorl, htrl, delu, delrho, delp, nx float64) (res float64) {
		var (
			delta            = arl / 20
			phi1, phi2, phi3 = phi(url-arl, delta), phi(url, delta), phi(url+arl, delta)
			ooarl2           = 1 / (arl * arl)
			f1               = (delp - rhorl*arl*delu) * 0.5 * ooarl2
			f2               = delrho - delp*ooarl2
			f3               = (delp + rhorl*arl*delu) * 0.5 * ooarl2
		)
		res = nx * (phi1*f1 + phi2*f2 + phi3*f3)
		return
	})
	eRhoU := URL.Copy().Apply8(aRL, RhoRL, HtRL, DelU, DelRho, DelP, el.NX, func(url, arl, rhorl, htrl, delu, delrho, delp, nx float64) (res float64) {
		var (
			delta            = arl / 20
			phi1, phi2, phi3 = phi(url-arl, delta), phi(url, delta), phi(url+arl, delta)
			ooarl2           = 1 / (arl * arl)
			f1               = (delp - rhorl*arl*delu) * 0.5 * ooarl2
			f2               = delrho - delp*ooarl2
			f3               = (delp + rhorl*arl*delu) * 0.5 * ooarl2
		)
		res = nx * (phi1*f1*(url-arl) + phi2*(f2*url+rhorl*delu) + phi3*f3*(url+arl))
		return
	})
	eEner := URL.Copy().Apply8(aRL, RhoRL, HtRL, DelU, DelRho, DelP, el.NX, func(url, arl, rhorl, htrl, delu, delrho, delp, nx float64) (res float64) {
		var (
			delta            = arl / 20
			phi1, phi2, phi3 = phi(url-arl, delta), phi(url, delta), phi(url+arl, delta)
			ooarl2           = 1 / (arl * arl)
			f1               = (delp - rhorl*arl*delu) * 0.5 * ooarl2
			f2               = delrho - delp*ooarl2
			f3               = (delp + rhorl*arl*delu) * 0.5 * ooarl2
		)
		res = nx * (phi1*f1*(htrl-arl*url) + phi2*(f2*url*url*0.5+rhorl*url*delu) + phi3*f3*(htrl+url*arl))
		return
	})

	fRho = fAve(RhoF).Subtract(eRho)
	fRhoU = fAve(RhoUF).Subtract(eRhoU)
	fEner = fAve(EnerF).Subtract(eEner)
	return
}

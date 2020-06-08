package Euler1D

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

type Euler struct {
	// Input parameters
	CFL, FinalTime  float64
	El, El_S        *DG1D.Elements1D
	RHSOnce         sync.Once
	State           *FieldState
	Rho, RhoU, Ener utils.Matrix
	In, Out         *State
	plotOnce        sync.Once
	chart           *chart2d.Chart2D
	colorMap        *utils2.ColorMap
	model           ModelType
	bc              BC_TYPE
	Case            CaseType
	frameCount      int
	useLimiter      bool
	FluxRanger      utils.R2
	FluxSubset      utils.Index
}

type CaseType uint

const (
	SOD_TUBE CaseType = iota
	DENSITY_WAVE
	COLLISION
	FREESTREAM
)

type ModelType uint

const (
	Galerkin_LF ModelType = iota
	DFR_LaxFriedrichs
	DFR_Roe
	DFR_Average
)

var (
	model_names = []string{
		"Galerkin Integration, Lax Flux",
		"DFR Integration, Lax Friedrichs Flux",
		"DFR Integration, Roe Flux",
		"DFR Integration, Average Flux",
	}
)

func NewEuler(CFL, FinalTime, XMax float64, N, K int, model ModelType, Case CaseType) (c *Euler) {
	switch Case {
	case DENSITY_WAVE:
		XMax = math.Max(XMax, 2)
	}
	VX, EToV := DG1D.SimpleMesh1D(0, XMax, K)
	c = &Euler{
		CFL:       CFL,
		State:     NewFieldState(),
		FinalTime: FinalTime,
		model:     model,
		Case:      Case,
	}
	switch model {
	case DFR_Roe, DFR_LaxFriedrichs, DFR_Average:
		c.El = DG1D.NewElements1D(N+2, VX, EToV)
		c.El_S = DG1D.NewElements1D(N, VX, EToV, DG1D.GAUSS)
	case Galerkin_LF:
		c.El = DG1D.NewElements1D(N, VX, EToV)
		c.El_S = c.El
	}
	c.MapSolutionSubset()
	c.State.Gamma = 1.4
	fmt.Printf("Euler Equations in 1 Dimension\n")
	switch c.Case {
	case DENSITY_WAVE:
		c.InitializeDWave()
		c.bc = PERIODIC
		fmt.Printf("Solving Density Wave\n")
	case SOD_TUBE:
		c.InitializeSOD()
		c.bc = RIEMANN
		c.useLimiter = true
		fmt.Printf("Solving Sod's Shock Tube\n")
	case FREESTREAM:
		c.InitializeFS()
		c.bc = RIEMANN
		fmt.Printf("Solving Freestream\n")
	case COLLISION:
		fallthrough
	default:
		fmt.Printf("Solving Shock Collision\n")
		c.InitializeSOD()
		c.bc = RIEMANN
		c.Out = c.In
		c.useLimiter = true
	}
	fmt.Printf("Algorithm: %s\n", model_names[c.model])
	if c.useLimiter {
		fmt.Printf("Solution is limited using SlopeLimit\n")
	} else {
		fmt.Printf("Solution is not limited\n")
	}
	fmt.Printf("CFL = %8.4f, Polynomial Degree N = %d (1 is linear), Num Elements K = %d\n\n\n", CFL, N, K)
	return
}

func (c *Euler) InitializeFS() {
	c.In = NewStateP(c.State.Gamma, 1, 0, 1)
	c.Out = c.In
	var (
		el = c.El_S
		FS = c.In
	)
	c.Rho = utils.NewMatrix(el.Np, el.K).AddScalar(FS.Rho)
	c.RhoU = utils.NewMatrix(el.Np, el.K).AddScalar(FS.RhoU)
	c.Ener = utils.NewMatrix(el.Np, el.K).AddScalar(FS.Ener)
}

func (c *Euler) InitializeSOD() {
	var (
		el                = c.El_S
		MassMatrix, VtInv utils.Matrix
		err               error
		npOnes            = utils.NewVectorConstant(el.Np, 1)
		s                 = c.State
	)
	c.In = NewStateP(c.State.Gamma, 1, 0, 1)
	c.Out = NewStateP(c.State.Gamma, 0.125, 0, 0.1)
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

func (c *Euler) InitializeDWave() {
	/*
		Rho(x,t) = 2 + sin(pi *(x-t))
		U = P = 1
	*/
	c.In = NewStateP(c.State.Gamma, 2, 2, 1)
	c.Out = c.In
	var (
		el = c.El_S
	)
	c.Rho = utils.NewMatrix(el.Np, el.K).Apply2(el.X, func(base, x float64) (rho float64) {
		rho = 2 + math.Sin(math.Pi*x)
		return
	})
	c.RhoU = c.Rho.Copy()
	// p = (e - q) * (fs.Gamma - 1)
	c.Ener = utils.NewMatrix(el.Np, el.K).Apply2(c.Rho, func(base, rho float64) (e float64) {
		e = 1./(c.State.Gamma-1.) + 0.5*rho
		return
	})
}

func (c *Euler) MapSolutionSubset() {
	/*
		In the DFR approach, there are Np solution points and Np+2 Flux points
		This index carves the solution points out of a full Np+2 space
	*/
	var (
		el  = c.El
		elS = c.El_S
	)
	c.FluxRanger = utils.NewR2(el.Np, el.K)
	if el.Np != elS.Np {
		c.FluxSubset = c.FluxRanger.Range("1:-1", ":")
	} else {
		c.FluxSubset = c.FluxRanger.Range(":", ":")
	}
	/*
		switch c.model {
		case DFR_LaxFriedrichs, DFR_Roe, DFR_Average:
			c.FluxSubset = c.FluxRanger.Range("1:-1", ":")
			//c.FluxSubset = c.FluxRanger.Range(":", ":")
		case Galerkin_LF:
			c.FluxSubset = c.FluxRanger.Range(":", ":")
		}
	*/
	return
}

func (c *Euler) Run(showGraph bool, graphDelay ...time.Duration) {
	var (
		el           = c.El
		elS          = c.El_S
		logFrequency = 50
		//s             = c.State
		rhs  func(rho, rhou, ener *utils.Matrix) (rhsRho, rhsRhoU, rhsEner utils.Matrix)
		iRho float64
	)
	xmin := elS.X.Row(1).Subtract(elS.X.Row(0)).Apply(math.Abs).Min()
	switch c.model {
	case Galerkin_LF:
		rhs = c.RHS_GK
	case DFR_Roe, DFR_LaxFriedrichs, DFR_Average:
		rhs = c.RHS_DFR
	}
	var Time, dt float64
	var tstep int
	for Time < c.FinalTime {
		/*
			Third Order Runge-Kutta time advancement
		*/
		// SSP RK Stage 1
		rhsRho, rhsRhoU, rhsEner := rhs(&c.Rho, &c.RhoU, &c.Ener)
		iRho = c.Plot(Time, showGraph, graphDelay)
		dt = c.CalculateDT(xmin, Time)
		update1 := func(u0, rhs float64) (u1 float64) {
			u1 = u0 + dt*rhs
			return
		}
		rho1 := c.Rho.Copy().Apply2(rhsRho, update1)
		rhou1 := c.RhoU.Copy().Apply2(rhsRhoU, update1)
		ener1 := c.Ener.Copy().Apply2(rhsEner, update1)

		// SSP RK Stage 2
		rhsRho, rhsRhoU, rhsEner = rhs(&rho1, &rhou1, &ener1)
		update2 := func(u0, u1, rhs float64) (u2 float64) {
			u2 = (3*u0 + u1 + rhs*dt) * (1. / 4.)
			return
		}
		rho2 := c.Rho.Copy().Apply3(rho1, rhsRho, update2)
		rhou2 := c.RhoU.Copy().Apply3(rhou1, rhsRhoU, update2)
		ener2 := c.Ener.Copy().Apply3(ener1, rhsEner, update2)

		// SSP RK Stage 3
		rhsRho, rhsRhoU, rhsEner = rhs(&rho2, &rhou2, &ener2)
		update3 := func(u0, u2, rhs float64) (u3 float64) {
			u3 = (u0 + 2*u2 + 2*dt*rhs) * (1. / 3.)
			return
		}
		c.Rho.Apply3(rho2, rhsRho, update3)
		c.RhoU.Apply3(rhou2, rhsRhoU, update3)
		c.Ener.Apply3(ener2, rhsEner, update3)

		Time += dt
		tstep++
		isDone := math.Abs(Time-c.FinalTime) < 0.000001
		if tstep%logFrequency == 0 || isDone {
			fmt.Printf("Time = %8.4f, max_resid[%d] = %8.4f, emin = %8.6f, emax = %8.6f\n", Time, tstep, rhsEner.Max(), c.Ener.Min(), c.Ener.Max())
			if isDone {
				switch c.Case {
				case SOD_TUBE:
					sod := sod_shock_tube.NewSOD(Time)
					x, rho, _, _, _ := sod.Get()
					fmt.Printf("SOD Shock Location = %5.4f\n", x[len(x)-2])
					iRho = integrate(x, rho)
					iRhoModel := integrate(elS.X.M.RawMatrix().Data, c.Rho.RawMatrix().Data)
					logErr := math.Log10(math.Abs(iRho - iRhoModel))
					rms_rho, rms_rhou, rms_e, max_rho, max_rhou, max_e := sod_error_calc(elS.X, c.Rho, c.RhoU, c.Ener, Time)
					if math.Abs(Time-c.FinalTime) < 0.001 {
						fmt.Printf("Rho Integration Check: Exact = %5.4f, Model = %5.4f, Log10 Error = %5.4f\n", iRho, iRhoModel, logErr)
						/*
							fmt.Printf("K=%v, N=%v, CFL=%v, RMS Errors: Log10 Rho, RhoU, E Error = %5.4f, %5.4f, %5.4f, MaxErr: Rho, RhoU, E = %5.4f, %5.4f, %5.4f\n",
								el.K, el.Np-1, c.CFL, math.Log10(rms_rho), math.Log10(rms_rhou), math.Log10(rms_e),
								math.Log10(max_rho), math.Log10(max_rhou), math.Log10(max_e))
						*/
						fmt.Printf("%s\n", "case,K,N,CFL,Log10_Rho_rms,Log10_Rhou_rms,Log10_e_rms,Log10_rho_max,Log10_rhou_max,Log10_e_max")
						fmt.Printf("\"%s\",%d,%d,%5.4f,%5.4f,%5.4f,%5.4f,%5.4f,%5.4f,%5.4f\n",
							model_names[c.model], el.K, el.Np-1, c.CFL, math.Log10(rms_rho), math.Log10(rms_rhou), math.Log10(rms_e),
							math.Log10(max_rho), math.Log10(max_rhou), math.Log10(max_e))
					}
				case DENSITY_WAVE:
					rms_rho, max_rho := dwaveErrorCalc(elS.X, c.Rho, Time)
					fmt.Printf("%s\n", "case,K,N,CFL,Log10_Rho_rms,Log10_rho_max")
					fmt.Printf("\"%s\",%d,%d,%5.4f,%5.4f,%5.4f\n",
						model_names[c.model], el.K, el.Np-1, c.CFL, math.Log10(rms_rho), math.Log10(max_rho))
				}
				if !showGraph {
					return
				}
			}
			if isDone && showGraph {
				for {
					time.Sleep(time.Second)
				}
			}
		}
	}
	return
}

func (c *Euler) CalculateDT(xmin, Time float64) (dt float64) {
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

type BC_TYPE uint

const (
	RIEMANN BC_TYPE = iota
	PERIODIC
)

func (c *Euler) RHS_DFR(Rhop, RhoUp, Enerp *utils.Matrix) (rhsRho, rhsRhoU, rhsEner utils.Matrix) {
	var (
		el                 = c.El
		elS                = c.El_S
		s                  = c.State
		fRho, fRhoU, fEner utils.Matrix
		Rho, RhoU, Ener    = *Rhop, *RhoUp, *Enerp
		RhoF, RhoUF, EnerF utils.Matrix
		limiter            = c.useLimiter
		slopeLimiterM      = 20.
	)
	if limiter {
		// Slope Limit the solution fields
		*Rhop = elS.SlopeLimitN(*Rhop, slopeLimiterM)
		*RhoUp = elS.SlopeLimitN(*RhoUp, slopeLimiterM)
		*Enerp = elS.SlopeLimitN(*Enerp, slopeLimiterM)
		Rho, RhoU, Ener = *Rhop, *RhoUp, *Enerp
	}
	RhoF, RhoUF, EnerF = s.Update(Rho, RhoU, Ener, c.FluxRanger, c.FluxSubset)
	if el.Np != elS.Np {
		c.CopyBoundary(RhoF)
		c.CopyBoundary(RhoUF)
		c.CopyBoundary(EnerF)
	}

	switch c.model {
	case DFR_Average:
		fRho, fRhoU, fEner = c.AveFlux(Rho, RhoU, Ener, RhoF, RhoUF, EnerF, elS.VmapM, elS.VmapP, el.VmapM, el.VmapP)
	case DFR_LaxFriedrichs:
		fRho, fRhoU, fEner = c.LaxFlux(Rho, RhoU, Ener, RhoF, RhoUF, EnerF, elS.VmapM, elS.VmapP, el.VmapM, el.VmapP)
	case DFR_Roe:
		fRho, fRhoU, fEner = c.RoeFlux(Rho, RhoU, Ener, RhoF, RhoUF, EnerF, elS.VmapM, elS.VmapP, el.VmapM, el.VmapP)
	}

	switch c.bc {
	case RIEMANN:
		c.RiemannBC_DFR(Rho, RhoU, Ener, RhoF, RhoUF, EnerF, &fRho, &fRhoU, &fEner)
	case PERIODIC:
		c.PeriodicBC_DFR(Rho, RhoU, Ener, RhoF, RhoUF, EnerF, elS.VmapI, elS.VmapO, el.VmapI, el.VmapO, &fRho, &fRhoU, &fEner)
	}

	//fmt.Println(RhoUF.Print("RhoUF Before Assign"))
	// Set face flux within global flux
	RhoF.AssignVector(el.VmapM, fRho)
	RhoUF.AssignVector(el.VmapM, fRhoU)
	EnerF.AssignVector(el.VmapM, fEner)

	/*
		fmt.Println(RhoUF.Print("RhoUF After Assign"))
		fmt.Println(el.Dr.Mul(RhoUF).Print("Dr*RhoUF"))
		fmt.Println(el.Dr.Mul(RhoUF).ElMul(el.Rx).Print("Dr*RhoUF.*Rx"))
		os.Exit(1)
	*/
	// Calculate RHS
	rhsRho = el.Dr.Mul(RhoF).Subset(c.FluxSubset, elS.Np, el.K).ElMul(elS.Rx).Scale(-1)
	rhsRhoU = el.Dr.Mul(RhoUF).Subset(c.FluxSubset, elS.Np, el.K).ElMul(elS.Rx).Scale(-1)
	rhsEner = el.Dr.Mul(EnerF).Subset(c.FluxSubset, elS.Np, el.K).ElMul(elS.Rx).Scale(-1)
	return
}

func (c *Euler) RHS_GK(Rhop, RhoUp, Enerp *utils.Matrix) (rhsRho, rhsRhoU, rhsEner utils.Matrix) {
	var (
		el                                                 = c.El
		nrF, ncF                                           = el.Nfp * el.NFaces, el.K
		s                                                  = c.State
		dRho, dRhoU, dEner, dRhoF, dRhoUF, dEnerF, LFcDiv2 utils.Matrix
		Rho, RhoU, Ener                                    = *Rhop, *RhoUp, *Enerp
		RhoF, RhoUF, EnerF                                 utils.Matrix
		limiter                                            = c.useLimiter
		slopeLimiterM                                      = 20.
	)
	if limiter {
		// Slope Limit the solution fields
		*Rhop = el.SlopeLimitN(*Rhop, slopeLimiterM)
		*RhoUp = el.SlopeLimitN(*RhoUp, slopeLimiterM)
		*Enerp = el.SlopeLimitN(*Enerp, slopeLimiterM)
		Rho, RhoU, Ener = *Rhop, *RhoUp, *Enerp
	}
	RhoF, RhoUF, EnerF = s.Update(Rho, RhoU, Ener, c.FluxRanger, c.FluxSubset)

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

	switch c.bc {
	case RIEMANN:
		c.RiemannBC(Rho, RhoU, Ener, RhoF, RhoUF, EnerF, &dRhoF, &dRhoUF, &dEnerF)
	default:
		panic("not implemented")
	}

	// RHS Computation
	rhsRho = el.LIFT.Mul(dRhoF.ElMul(el.FScale)).Subtract(el.Dr.Mul(RhoF).ElMul(el.Rx))
	rhsRhoU = el.LIFT.Mul(dRhoUF.ElMul(el.FScale)).Subtract(el.Dr.Mul(RhoUF).ElMul(el.Rx))
	rhsEner = el.LIFT.Mul(dEnerF.ElMul(el.FScale)).Subtract(el.Dr.Mul(EnerF).ElMul(el.Rx))
	return
}

func (c *Euler) PeriodicBC_DFR(Rho, RhoU, Ener, RhoF, RhoUF, EnerF utils.Matrix, vmapIS, vmapOS, vmapI, vmapO utils.Index, dRhoF, dRhoUF, dEnerF *utils.Matrix) {
	// Periodic Boundary condition
	fRho, fRhoU, fEner := c.RoeFlux(Rho, RhoU, Ener, RhoF, RhoUF, EnerF, vmapIS, vmapOS, vmapI, vmapO)
	/*
		fmt.Println(dRhoF.Print("dRhoF before"))
		fmt.Println(dRhoUF.Print("dRhoUF before"))
		fmt.Println(dEnerF.Print("dEnerF before"))
	*/
	dRhoF.AssignVector(c.El.MapI, fRho)
	dRhoUF.AssignVector(c.El.MapI, fRhoU)
	dEnerF.AssignVector(c.El.MapI, fEner)
	dRhoF.AssignVector(c.El.MapO, fRho)
	dRhoUF.AssignVector(c.El.MapO, fRhoU)
	dEnerF.AssignVector(c.El.MapO, fEner)
	/*
		fmt.Println(dRhoF.Print("dRhoF after"))
		fmt.Println(dRhoUF.Print("dRhoUF after"))
		fmt.Println(dEnerF.Print("dEnerF after"))
		os.Exit(1)
	*/
}

func (c *Euler) RiemannBC(Rho, RhoU, Ener, RhoF, RhoUF, EnerF utils.Matrix, dRhoF, dRhoUF, dEnerF *utils.Matrix) {
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

func (c *Euler) RiemannBC_DFR(Rho, RhoU, Ener, RhoF, RhoUF, EnerF utils.Matrix, dRhoF, dRhoUF, dEnerF *utils.Matrix) {
	var (
		s   = c.State
		el  = c.El
		elS = c.El_S
		// Sod's problem: Shock tube with jump in middle
		In  = c.In
		Out = c.Out
	)
	// Boundary conditions for Sod's problem
	// Inflow
	lmI := s.LM.SubsetVector(elS.VmapI).Scale(0.5)
	nxI := el.NX.SubsetVector(el.MapI)
	bFunc_dfr(dRhoF, Rho, RhoF, lmI, nxI, In.Rho, In.RhoF, el.MapI, elS.VmapI, el.VmapI)
	bFunc_dfr(dRhoUF, RhoU, RhoUF, lmI, nxI, In.RhoU, In.RhoUF, el.MapI, elS.VmapI, el.VmapI)
	bFunc_dfr(dEnerF, Ener, EnerF, lmI, nxI, In.Ener, In.EnerF, el.MapI, elS.VmapI, el.VmapI)
	// Outflow
	lmO := s.LM.SubsetVector(elS.VmapO).Scale(0.5)
	nxO := el.NX.SubsetVector(el.MapO)
	bFunc_dfr(dRhoF, Rho, RhoF, lmO, nxO, Out.Rho, Out.RhoF, el.MapO, elS.VmapO, el.VmapO)
	bFunc_dfr(dRhoUF, RhoU, RhoUF, lmO, nxO, Out.RhoU, Out.RhoUF, el.MapO, elS.VmapO, el.VmapO)
	bFunc_dfr(dEnerF, Ener, EnerF, lmO, nxO, Out.Ener, Out.EnerF, el.MapO, elS.VmapO, el.VmapO)
}

func (c *Euler) Plot(timeT float64, showGraph bool, graphDelay []time.Duration) (iRho float64) {
	var (
		el         = c.El
		elS        = c.El_S
		fmin, fmax float32
	)
	switch c.Case {
	case SOD_TUBE, COLLISION:
		fmin, fmax = float32(-0.1), float32(2.6)
	case DENSITY_WAVE:
		fmin, fmax = float32(1.0), float32(4.0)
	case FREESTREAM:
		fmin, fmax = float32(-0.1), float32(2.6)
	}
	if !showGraph {
		return
	}
	c.plotOnce.Do(func() {
		c.chart = chart2d.NewChart2D(1920, 1280, float32(el.X.Min()), float32(el.X.Max()), fmin, fmax)
		c.colorMap = utils2.NewColorMap(-1, 1, 1)
		go c.chart.Plot()
	})
	pSeries := func(field utils.Matrix, name string, color float32, gl chart2d.GlyphType) {
		var (
			x utils.Matrix
		)
		x = el.X
		if elS.Np != el.Np {
			x = elS.X
		}
		if err := c.chart.AddSeries(name, x.Transpose().RawMatrix().Data, field.Transpose().RawMatrix().Data,
			gl, chart2d.Solid, c.colorMap.GetRGB(color)); err != nil {
			panic("unable to add graph series")
		}
	}
	pSeries(c.Rho, "Rho", -0.7, chart2d.NoGlyph)
	pSeries(c.RhoU, "RhoU", 0.0, chart2d.NoGlyph)
	pSeries(c.Ener, "Ener", 0.7, chart2d.NoGlyph)
	c.frameCount++
	check := int(math.Log10(float64(el.K * el.Np / 5)))
	if c.frameCount%check == 0 || math.Abs(timeT-c.FinalTime) < 0.001 {
		switch c.Case {
		case SOD_TUBE:
			iRho = AddAnalyticSod(c.chart, c.colorMap, timeT)
		case DENSITY_WAVE:
			AddAnalyticDWave(c.chart, c.colorMap, elS.X, timeT)
		}
	}
	if len(graphDelay) != 0 {
		time.Sleep(graphDelay[0])
	}
	return
}

func DWaveCalc(X utils.Matrix, timeT float64) (x, rho []float64) {
	Rho := X.Copy().Apply(func(x float64) (rho float64) {
		rho = 2 + math.Sin(math.Pi*(x-timeT))
		return
	})
	x, rho = X.RawMatrix().Data, Rho.RawMatrix().Data
	return
}

func AddAnalyticDWave(chart *chart2d.Chart2D, colorMap *utils2.ColorMap, X utils.Matrix, timeT float64) {
	x, rho := DWaveCalc(X, timeT)
	if err := chart.AddSeries("ExactRho", x, rho, chart2d.XGlyph, chart2d.NoLine, colorMap.GetRGB(-0.7)); err != nil {
		panic("unable to add exact solution Rho")
	}
	return
}

func dwaveErrorCalc(X, Rho utils.Matrix, t float64) (rms_rho, max_rho float64) {
	var (
		RhoData = Rho.RawMatrix().Data
	)
	_, rhoD := DWaveCalc(X, t)
	for i, rho := range RhoData {
		rho_err := utils.POW(rho-rhoD[i], 2)
		rms_rho += rho_err
		rho_err = math.Sqrt(rho_err)
		max_rho = math.Max(rho_err, max_rho)
	}
	rms_rho = math.Sqrt(rms_rho / float64(len(RhoData)))
	return
}

func AddAnalyticSod(chart *chart2d.Chart2D, colorMap *utils2.ColorMap, timeT float64) (iRho float64) {
	sod := sod_shock_tube.NewSOD(timeT)
	X, Rho, _, RhoU, E := sod.Get()
	if err := chart.AddSeries("ExactRho", X, Rho, chart2d.XGlyph, chart2d.NoLine, colorMap.GetRGB(-0.7)); err != nil {
		panic("unable to add exact solution Rho")
	}
	if err := chart.AddSeries("ExactRhoU", X, RhoU, chart2d.XGlyph, chart2d.NoLine, colorMap.GetRGB(0.0)); err != nil {
		panic("unable to add exact solution RhoU")
	}
	if err := chart.AddSeries("ExactE", X, E, chart2d.XGlyph, chart2d.NoLine, colorMap.GetRGB(0.7)); err != nil {
		panic("unable to add exact solution E")
	}
	iRho = integrate(X, Rho)
	return
}

func integrate(x, u []float64) (result float64) {
	L := len(x)
	for i := 0; i < L-1; i++ {
		delx := x[i+1] - x[i]
		uave := 0.5 * (u[i+1] + u[i])
		result += uave * delx
	}
	return
}

func sod_error_calc(X, Rho, RhoU, E utils.Matrix, t float64) (rms_rho, rms_rhou, rms_e, max_rho, max_rhou, max_e float64) {
	var (
		Xdata    = X.RawMatrix().Data
		RhoData  = Rho.RawMatrix().Data
		RhoUData = RhoU.RawMatrix().Data
		EData    = E.RawMatrix().Data
	)
	sod := sod_shock_tube.NewSOD(t)
	for i, x := range Xdata {
		// Only validate using flow left of center, excluding the shock and contact discontinuity
		if x < 0.5 && x > 0.05 {
			sod_rho, _, _, sod_e, sod_rhou := sod.Getx(x)
			rho, rhou, e := RhoData[i], RhoUData[i], EData[i]
			rho_err := utils.POW(rho-sod_rho, 2)
			rhou_err := utils.POW(rhou-sod_rhou, 2)
			e_err := utils.POW(e-sod_e, 2)
			rms_rho += rho_err
			rms_rhou += rhou_err
			rms_e += e_err
			rho_err, rhou_err, e_err = math.Sqrt(rho_err), math.Sqrt(rhou_err), math.Sqrt(e_err)
			max_rho = math.Max(rho_err, max_rho)
			max_rhou = math.Max(rhou_err, max_rhou)
			max_e = math.Max(e_err, max_e)
		}
	}
	rms_rho = math.Sqrt(rms_rho / float64(len(Xdata)))
	rms_rhou = math.Sqrt(rms_rhou / float64(len(Xdata)))
	rms_e = math.Sqrt(rms_e / float64(len(Xdata)))
	return
}

func (c *Euler) AveFlux(Rho, RhoU, Ener, RhoF, RhoUF, EnerF utils.Matrix, vmapMS, vmapPS, vmapM, vmapP utils.Index) (fRho, fRhoU, fEner utils.Matrix) {
	var (
		el       = c.El
		nrF, ncF = el.Nfp * el.NFaces, el.K
	)
	// Compute Lax-Friedrichs flux
	// Face flux average
	fAve := func(U utils.Matrix) (Uavg utils.Matrix) {
		Uavg = U.Subset(vmapM, nrF, ncF).Add(U.Subset(vmapP, nrF, ncF)).Scale(0.5)
		return
	}
	// Compute numerical flux at faces
	fRho = fAve(RhoF)
	fRhoU = fAve(RhoUF)
	fEner = fAve(EnerF)
	return
}

func (c *Euler) LaxFlux(Rho, RhoU, Ener, RhoF, RhoUF, EnerF utils.Matrix, vmapMS, vmapPS, vmapM, vmapP utils.Index) (fRho, fRhoU, fEner utils.Matrix) {
	var (
		el       = c.El
		s        = c.State
		nrF, ncF = el.Nfp * el.NFaces, el.K
	)
	// Compute Lax-Friedrichs flux
	// Face jumps in primary and flux variables
	fJump := func(U utils.Matrix) (dU utils.Matrix) {
		dU = U.Subset(vmapMS, nrF, ncF).Subtract(U.Subset(vmapPS, nrF, ncF)).ElMul(el.NX)
		return
	}
	// Face flux average
	fAve := func(U utils.Matrix) (Uavg utils.Matrix) {
		Uavg = U.Subset(vmapM, nrF, ncF).Add(U.Subset(vmapP, nrF, ncF)).Scale(0.5)
		return
	}
	// Max eigenvalue
	LFc := s.LM.Subset(vmapMS, nrF, ncF).Apply2(s.LM.Subset(vmapPS, nrF, ncF), math.Max)
	// Compute numerical flux at faces
	fRho = fAve(RhoF).Add(fJump(Rho).ElMul(LFc).Scale(0.5))
	fRhoU = fAve(RhoUF).Add(fJump(RhoU).ElMul(LFc).Scale(0.5))
	fEner = fAve(EnerF).Add(fJump(Ener).ElMul(LFc).Scale(0.5))
	return
}

func (c *Euler) CopyBoundary(U utils.Matrix) {
	var (
		el = c.El
	)
	/*
		Copy the flux values from the interior to the edge in prep for construction
	*/
	U.M.SetRow(0, U.Row(1).RawVector().Data)
	U.M.SetRow(el.Np-1, U.Row(el.Np-2).RawVector().Data)
}
func (c *Euler) RoeFlux(Rho, RhoU, Ener, RhoF, RhoUF, EnerF utils.Matrix, vmapMS, vmapPS, vmapM, vmapP utils.Index) (fRho, fRhoU, fEner utils.Matrix) {
	var (
		el       = c.El
		s        = c.State
		nrF, ncF int
	)
	nrF, ncF = el.Nfp*el.NFaces, len(vmapM)/2
	if ncF == 0 {
		ncF = 1
	}
	fL := func(U utils.Matrix) (Ul utils.Matrix) {
		Ul = U.Subset(vmapMS, nrF, ncF)
		return
	}
	fR := func(U utils.Matrix) (Ul utils.Matrix) {
		Ul = U.Subset(vmapPS, nrF, ncF)
		return
	}
	// Face jumps in primary and flux variables
	fJump := func(U utils.Matrix) (dU utils.Matrix) {
		dU = U.Subset(vmapPS, nrF, ncF).Subtract(U.Subset(vmapMS, nrF, ncF))
		return
	}
	// Face average
	fAve := func(U utils.Matrix) (Uavg utils.Matrix) {
		Uavg = U.Subset(vmapP, nrF, ncF).Add(U.Subset(vmapM, nrF, ncF)).Scale(0.5)
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
		res = nx * (phi1*f1*(url-arl) + phi2*(f2*url) + phi3*f3*(url+arl))
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
		res = nx * (phi1*f1*(htrl-arl*url) + phi2*(f2*url*url*0.5) + phi3*f3*(htrl+url*arl))
		return
	})

	fRho = fAve(RhoF).Subtract(eRho)
	fRhoU = fAve(RhoUF).Subtract(eRhoU)
	fEner = fAve(EnerF).Subtract(eEner)
	return
}

type FieldState struct {
	Gamma                    float64
	U, Q, Pres, CVel, LM, Ht utils.Matrix
	Temp                     utils.Matrix
}

func NewFieldState() (fs *FieldState) {
	return &FieldState{}
}

func (fs *FieldState) Update(Rho, RhoU, Ener utils.Matrix, FluxRanger utils.R2, FluxSubset utils.Index) (RhoF, RhoUF, EnerF utils.Matrix) {
	/*
		Calorically Perfect Gas, R = 1
	*/
	var (
		Gamma = 1.4
		Cv    = 1. / (Gamma - 1.)
		Cp    = Gamma * Cv
	)
	fs.U = RhoU.Copy().ElDiv(Rho) // Velocity
	fs.Q = fs.U.Copy().Apply2(Rho, func(u, rho float64) (q float64) {
		q = 0.5 * u * u * rho
		return
	})
	// (gamma-1.0)*(Ener - 0.5*(rhou).^2./rho)
	fs.Pres = Ener.Copy().Apply2(fs.Q, func(e, q float64) (p float64) {
		p = (e - q) * (fs.Gamma - 1)
		return
	})
	// sqrt(gamma*pres./rho)
	fs.CVel = fs.Pres.Copy().Apply2(Rho, func(p, rho float64) (cvel float64) {
		cvel = math.Sqrt(fs.Gamma * p / rho)
		return
	})
	//  abs(rhou./rho)+cvel
	fs.LM = fs.U.Copy().Apply2(fs.CVel, func(u, cvel float64) (lm float64) {
		lm = math.Abs(u) + cvel
		return
	})
	//	Temp = (Ener - 0.5*(rhou).^2./rho)./rho
	fs.Temp = Ener.Copy().Apply3(fs.Q, Rho, func(e, q, rho float64) (temp float64) {
		temp = (e - q) / rho
		return
	})
	fs.Ht = fs.Q.Copy().Apply3(Rho, fs.Temp, func(q, rho, temp float64) (ht float64) {
		ht = Cp * (q/rho + temp)
		return
	})
	fs.U.SetReadOnly("U")
	fs.Q.SetReadOnly("Q")
	fs.Pres.SetReadOnly("Pres")
	fs.CVel.SetReadOnly("CVel")
	fs.LM.SetReadOnly("LM")
	fs.Ht.SetReadOnly("Ht")
	RhoF = utils.NewMatrix(FluxRanger.Dims())
	RhoUF = utils.NewMatrix(FluxRanger.Dims())
	EnerF = utils.NewMatrix(FluxRanger.Dims())
	RhoF.AssignVector(FluxSubset, RhoU)
	RhoUF.AssignVector(FluxSubset, fs.Q.Copy().Apply2(fs.Pres, func(q, pres float64) (rhouf float64) {
		rhouf = 2*q + pres
		return
	}))
	EnerF.AssignVector(FluxSubset, Ener.Copy().Apply3(fs.Pres, fs.U, func(ener, pres, u float64) (enerf float64) {
		enerf = u * (ener + pres)
		return
	}))
	return
}

func (fs *FieldState) Print() {
	fmt.Println(fs.U.Print("U"))
	fmt.Println(fs.Q.Print("Q"))
	fmt.Println(fs.Pres.Print("Pres"))
	fmt.Println(fs.CVel.Print("CVel"))
	fmt.Println(fs.LM.Print("LM"))
	fmt.Println(fs.Ht.Print("Ht"))
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

func bFunc_dfr(dUF *utils.Matrix, U, UF utils.Matrix, lm, nx utils.Vector, uIO, ufIO float64, mapi, vmaps, vmap utils.Index) {
	// Characteristic BC using freestream conditions at boundary
	var (
		dUFVec utils.Matrix
	)
	dUFVec = nx.Outer(UF.SubsetVector(vmap).Subtract(utils.NewVectorConstant(len(vmap), ufIO))).Scale(0.5)
	dUFVec.Subtract(lm.Outer(U.SubsetVector(vmaps).Subtract(utils.NewVectorConstant(len(vmaps), uIO))))
	//dUF.AssignVector(mapi, dUFVec)
	nxI := nx.AtVec(0)
	dUF.AssignVector(mapi, dUFVec.Scale(-nxI*0.5).Add(UF.Subset(vmap, 1, 1)))
	return
}

type State struct {
	Gamma, Rho, RhoU, Ener float64
	RhoF, RhoUF, EnerF     float64
}

func NewState(gamma, rho, rhoU, ener float64) (s *State) {
	q := 0.5 * utils.POW(rhoU, 2) / rho
	p := (ener - q) * (gamma - 1.)
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
	q := 0.5 * rhoU * rhoU / rho
	ener := (p - q) / (gamma - 1.)
	return NewState(gamma, rho, rhoU, ener)
}

func (s *State) Print() (o string) {
	o = fmt.Sprintf("Rho = %v\nP = %v\nE = %v\nRhoU = %v\nRhoUF = %v\n",
		s.Rho, s.Ener*(s.Gamma-1.), s.Ener, s.RhoU, s.RhoUF)
	return
}

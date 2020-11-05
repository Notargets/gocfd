package Euler2D

import (
	"fmt"
	"math"

	"github.com/notargets/gocfd/DG2D"

	"github.com/notargets/avs/chart2d"
	utils2 "github.com/notargets/avs/utils"
	"github.com/notargets/gocfd/utils"
)

/*
	In the DFR scheme, we have two sets of points:
		- Solution points (inside element)
		- Flux Points (on element faces)
	The flux element operates on both sets of points, the solution element only operates on the internal solution points
*/
type Euler struct {
	// Input parameters
	MeshFile              string
	CFL, Gamma, FinalTime float64
	dfr                   *DG2D.DFR2D
	Q                     [4]utils.Matrix // Solution variables, stored at solution point locations, Np_solution x K
	Q_Face                [4]utils.Matrix // Solution variables, interpolated to and stored at edge point locations, Np_edge x K
	F_RT_DOF              [4]utils.Matrix // Normal Projected Flux, stored at flux/solution point locations, Np_flux x K
	chart                 *chart2d.Chart2D
	colorMap              *utils2.ColorMap
	model                 ModelType
	Case                  CaseType
}

type CaseType uint

const (
	FREESTREAM CaseType = iota
)

type ModelType uint

const (
	FLUX_LaxFriedrichs ModelType = iota
	FLUX_Roe
	FLUX_Average
)

var (
	modelNames = []string{
		"Lax Friedrichs Flux",
		"Roe Flux",
		"Average Flux",
	}
)

func NewEuler(CFL, FinalTime float64, N int, meshFile string, model ModelType, Case CaseType, verbose bool) (c *Euler) {
	c = &Euler{
		MeshFile:  meshFile,
		CFL:       CFL,
		FinalTime: FinalTime,
		model:     model,
		Case:      Case,
		Gamma:     1.4,
	}
	c.dfr = DG2D.NewDFR2D(N, meshFile)
	c.InitializeMemory()
	if verbose {
		fmt.Printf("Euler Equations in 2 Dimensions\n")
	}
	switch c.Case {
	case FREESTREAM:
		c.InitializeFS()
		if verbose {
			fmt.Printf("Solving Freestream\n")
		}
	default:
		panic("unknown case type")
	}
	if verbose {
		fmt.Printf("Calling avg flux for file: %s ...\n", c.MeshFile)
		fmt.Printf("done\n")
		fmt.Printf("Algorithm: %s\n", modelNames[c.model])
		fmt.Printf("CFL = %8.4f, Polynomial Degree N = %d (1 is linear), Num Elements K = %d\n\n\n", CFL, N, c.dfr.K)
	}
	return
}

func (c *Euler) AssembleRTNormalFlux() {
	/*
		Solver approach:
		0) Solution is stored on sol points as Q
		0a) Flux is computed and stored in X, Y component projections in the 2*Nint front of F_RT_DOF
		1) Solution is extrapolated to edge points in Q_Face from Q
		2) Edges are traversed, flux is calculated and projected onto edge face normals, scaled and placed into F_RT_DOF
	*/
	Kmax := c.dfr.K
	Nint := c.dfr.FluxElement.Nint
	// Calculate flux and project into R and S (transformed) directions for the internal points
	for i := 0; i < Nint; i++ {
		for k := 0; k < c.dfr.K; k++ {
			ind := k + i*Kmax
			Fr, Fs := c.CalculateFluxTransformed(k, i, c.Q)
			for n := 0; n < 4; n++ {
				rtD := c.F_RT_DOF[n].Data()
				rtD[ind], rtD[ind+Nint*Kmax] = Fr[n], Fs[n]
			}
		}
	}
	// Interpolate from solution points to edges using precomputed interpolation matrix
	for n := 0; n < 4; n++ {
		c.Q_Face[n] = c.dfr.FluxInterpMatrix.Mul(c.Q[n])
	}
	c.SetNormalFluxOnEdges()
	// TODO: Implement calculation of flux for BCs
}

func (c *Euler) SetNormalFluxOnEdges() {
	var (
		dfr   = c.dfr
		Nedge = dfr.FluxElement.Nedge
		Kmax  = dfr.K
	)
	norm := func(vec [2]float64) (n float64) {
		n = math.Sqrt(vec[0]*vec[0] + vec[1]*vec[1])
		return
	}
	normalize := func(vec [2]float64) (normed [2]float64) {
		n := norm(vec)
		for i := 0; i < 2; i++ {
			normed[i] = vec[i] / n
		}
		return
	}
	getScaledNormal := func(conn int, e *DG2D.Edge, en DG2D.EdgeNumber) (scaledNormal [2]float64) {
		revDir := bool(e.ConnectedTriDirection[conn])
		x1, x2 := dfr.Tris.GetEdgeCoordinates(en, revDir, dfr.VX, dfr.VY)
		dx := [2]float64{x2[0] - x1[0], x2[1] - x1[1]}
		scaledNormal = normalize([2]float64{-dx[1], dx[0]})
		scaledNormal[0] *= e.IInII[conn]
		scaledNormal[1] *= e.IInII[conn]
		return
	}
	// Handle only edges with two connected tris
	for en, e := range dfr.Tris.Edges {
		if e.NumConnectedTris == 2 {
			var (
				kL, kR                          = int(e.ConnectedTris[0]), int(e.ConnectedTris[1])
				edgeNumberL, edgeNumberR        = int(e.ConnectedTriEdgeNumber[0]), int(e.ConnectedTriEdgeNumber[1])
				shiftL, shiftR                  = edgeNumberL * Nedge, edgeNumberR * Nedge
				fluxLeft, fluxRight, fluxAve    [2][4]float64
				normalFluxLeft, normalFluxRight [4]float64
			)
			// Get scaling factor ||n|| for each edge, multiplied by untransformed normals
			normLeft, normRight := getScaledNormal(0, e, en), getScaledNormal(1, e, en)
			for i := 0; i < Nedge; i++ {
				iL := i + shiftL
				iR := Nedge - 1 - i + shiftR // Shared edges run in reverse order relative to each other
				fluxLeft[0], fluxLeft[1] = c.CalculateFlux(kL, iL, c.Q_Face)
				fluxRight[0], fluxRight[1] = c.CalculateFlux(kR, iR, c.Q_Face) // Reverse the right edge to match
				// Implement average flux (for now), later we will add Roe, Lax, etc flux calculations
				for n := 0; n < 4; n++ {
					for ii := 0; ii < 2; ii++ {
						fluxAve[ii][n] = 0.5 * (fluxLeft[ii][n] + fluxRight[ii][n])
					}
				}
				// Project the flux onto the scaled normal for each of left/right
				for n := 0; n < 4; n++ {
					normalFluxLeft[n] = normLeft[0]*fluxAve[0][n] + normLeft[1]*fluxAve[1][n]
					normalFluxRight[n] = normRight[0]*fluxAve[0][n] + normRight[1]*fluxAve[1][n]
				}
				// Place normed/scaled flux into the RT element space
				for n := 0; n < 4; n++ {
					rtD := c.F_RT_DOF[n].Data()
					indL := kL + iL*Kmax
					rtD[indL] = normalFluxLeft[n]
					indR := kR + iR*Kmax
					rtD[indR] = normalFluxLeft[n]
				}
			}
		}
	}
	return
}

/*
func (c *Euler) AverageFlux() {
	var (
		el    = c.dfr.FluxElement
		Nedge = el.Nedge
	)
	for _, e := range c.dfr.Tris.Edges {
		if e.BCType == DG2D.BC_None && e.NumConnectedTris == 2 {
			// We construct a shared flux
			k1, k2 := int(e.ConnectedTris[0]), int(e.ConnectedTris[1])
			for ii := 0; ii < 4; ii++ {
				fxD, fyD := c.Fx[ii].Data(), c.Fy[ii].Data()
				ind1, ind2 := c.EdgeStart(k1, e, 0), c.EdgeStart(k2, e, 1)
				for i := 0; i < Nedge; i++ {
					ind := i + ind1
					indR := (Nedge - 1) - i + ind2
					// Reverse the right relative to the left
					aveFx := 0.5 * (fxD[ind] + fxD[indR])
					aveFy := 0.5 * (fyD[ind] + fyD[indR])
					fxD[ind], fxD[indR] = aveFx, aveFx
					fyD[ind], fyD[indR] = aveFy, aveFy
				}
			}
		}
	}
}
*/

func (c *Euler) EdgeStart(k int, e *DG2D.Edge, conn int) (index int) {
	/*
			Flux points are stored as (KxNp) for each Flux
		    Within Np, the flux points are layed out like:
			<---- Nint ----><---- Nint ----><---Nedge----><---Nedge----><---Nedge---->
			         Solution Points          Edge 0 pts	Edge 1 pts	  Edge 2 pts
			<---- Nint ----><---- Nint ----><---Nedge----><---Nedge----><---Nedge---->
	*/
	var (
		el      = c.dfr.FluxElement
		Np      = el.Np
		Nint    = el.Nint
		Nedge   = el.Nedge
		edgeNum = e.ConnectedTriEdgeNumber[conn].Index()
	)
	index = k*Np + 2*Nint + edgeNum*Nedge
	return
}

func (c *Euler) InitializeFS() {
	var (
		rho, u, v, p = 1., 0., 0., 1. // Freestream state
		K            = c.dfr.K
		Np           = c.dfr.SolutionElement.Np
		Gamma        = c.Gamma
		GM1          = Gamma - 1 // R / Cv
	)
	q := 0.5 * rho * (u*u + v*v)
	rhoE := p/GM1 + q
	c.Q[0] = utils.NewMatrix(Np, K).AddScalar(rho)
	c.Q[1] = utils.NewMatrix(Np, K).AddScalar(rho * u)
	c.Q[2] = utils.NewMatrix(Np, K).AddScalar(rho * v)
	c.Q[3] = utils.NewMatrix(Np, K).AddScalar(rhoE)
}

func (c *Euler) InitializeMemory() {
	var (
		K      = c.dfr.K
		Nedge  = c.dfr.FluxElement.Nedge
		NpFlux = c.dfr.FluxElement.Np
	)
	for n := 0; n < 4; n++ {
		c.Q_Face[n] = utils.NewMatrix(Nedge*3, K)
		c.F_RT_DOF[n] = utils.NewMatrix(NpFlux, K)
	}
}

func (c *Euler) CalculateFluxTransformed(k, i int, Q [4]utils.Matrix) (Fr, Fs [4]float64) {
	var (
		J, _, Jdet = c.dfr.GetJacobian(k)
	)
	Fx, Fy := c.CalculateFlux(k, i, Q)
	for n := 0; n < 4; n++ {
		Fr[n] = (1. / Jdet) * (J[0]*Fx[n] + J[1]*Fy[n])
		Fs[n] = (1. / Jdet) * (J[2]*Fx[n] + J[3]*Fy[n])
	}
	return
}

func (c *Euler) CalculateFlux(k, i int, Q [4]utils.Matrix) (Fx, Fy [4]float64) {
	// From https://www.theoretical-physics.net/dev/fluid-dynamics/euler.html
	var (
		Gamma                 = c.Gamma
		GM1                   = Gamma - 1 // R / Cv
		q0D, q1D, q2D, q3D    = Q[0].Data(), Q[1].Data(), Q[2].Data(), Q[3].Data()
		Kmax                  = c.dfr.K
		ind                   = k + i*Kmax
		rho, rhoU, rhoV, rhoE = q0D[ind], q1D[ind], q2D[ind], q3D[ind]
	)
	u := rhoU / rho
	v := rhoV / rho
	u2 := u*u + v*v
	q := 0.5 * rho * u2
	p := GM1 * (rhoE - q)
	Fx, Fy =
		[4]float64{rhoU, rhoU*u + p, rhoU * v, u * (rhoE + p)},
		[4]float64{rhoV, rhoV * u, rhoV*v + p, v * (rhoE + p)}
	return
}

/*
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
		//	Third Order Runge-Kutta time advancement
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
						fmt.Printf("%s\n", "case,K,N,CFL,Log10_Rho_rms,Log10_Rhou_rms,Log10_e_rms,Log10_rho_max,Log10_rhou_max,Log10_e_max")
						fmt.Printf("\"%s\",%d,%d,%5.4f,%5.4f,%5.4f,%5.4f,%5.4f,%5.4f,%5.4f\n",
							modelNames[c.model], el.K, elS.Np-1, c.CFL, math.Log10(rms_rho), math.Log10(rms_rhou), math.Log10(rms_e),
							math.Log10(max_rho), math.Log10(max_rhou), math.Log10(max_e))
					}
				case DENSITY_WAVE:
					rms_rho, max_rho := dwaveErrorCalc(elS.X, c.Rho, Time)
					fmt.Printf("%s\n", "case,K,N,CFL,Log10_Rho_rms,Log10_rho_max")
					fmt.Printf("\"%s\",%d,%d,%5.4f,%5.4f,%5.4f\n",
						modelNames[c.model], el.K, el.Np-1, c.CFL, math.Log10(rms_rho), math.Log10(max_rho))
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

func (c *Euler) RHS_DFR(Rhop, RhoUp, Enerp *utils.Matrix) (rhsRho, rhsRhoU, rhsEner utils.Matrix) {
	var (
		el                          = c.El
		elS                         = c.El_S
		s                           = c.State
		fRho, fRhoU, fEner          utils.Matrix
		Rho, RhoU, Ener             = *Rhop, *RhoUp, *Enerp
		RhoF, RhoUF, EnerF          utils.Matrix
		limiter                     = c.useLimiter
		RhoFull, RhoUFull, EnerFull utils.Matrix
	)

	RhoFull, RhoUFull, EnerFull, RhoF, RhoUF, EnerF = s.Update(Rho, RhoU, Ener, c)
	if limiter {
		// Slope Limit the solution fields
		*Rhop = RhoFull.Subset(c.FluxSubset, elS.Np, elS.K)
		*RhoUp = RhoUFull.Subset(c.FluxSubset, elS.Np, elS.K)
		*Enerp = EnerFull.Subset(c.FluxSubset, elS.Np, elS.K)
	}

	switch c.model {
	case DFR_Average:
		fRho, fRhoU, fEner = c.AveFlux(RhoFull, RhoUFull, EnerFull, RhoF, RhoUF, EnerF, el.VmapM, el.VmapP)
	case DFR_LaxFriedrichs:
		fRho, fRhoU, fEner = c.LaxFlux(RhoFull, RhoUFull, EnerFull, RhoF, RhoUF, EnerF, el.VmapM, el.VmapP)
	case DFR_Roe:
		fRho, fRhoU, fEner = c.RoeFlux(RhoFull, RhoUFull, EnerFull, RhoF, RhoUF, EnerF, el.VmapM, el.VmapP)
	}

	switch c.bc {
	case RIEMANN:
		c.RiemannBC_DFR(RhoFull, RhoUFull, EnerFull, RhoF, RhoUF, EnerF, &fRho, &fRhoU, &fEner)
	case PERIODIC:
		c.PeriodicBC_DFR(RhoFull, RhoUFull, EnerFull, RhoF, RhoUF, EnerF, el.VmapI, el.VmapO, &fRho, &fRhoU, &fEner)
	}

	// Set face flux within global flux
	RhoF.AssignVector(el.VmapM, fRho)
	RhoUF.AssignVector(el.VmapM, fRhoU)
	EnerF.AssignVector(el.VmapM, fEner)

	// Calculate RHS
	rhsRho = el.Dr.Mul(RhoF).Subset(c.FluxSubset, elS.Np, el.K).ElMul(elS.Rx).Scale(-1)
	rhsRhoU = el.Dr.Mul(RhoUF).Subset(c.FluxSubset, elS.Np, el.K).ElMul(elS.Rx).Scale(-1)
	rhsEner = el.Dr.Mul(EnerF).Subset(c.FluxSubset, elS.Np, el.K).ElMul(elS.Rx).Scale(-1)
	return
}

func (c *Euler) PeriodicBC_DFR(Rho, RhoU, Ener, RhoF, RhoUF, EnerF utils.Matrix, vmapI, vmapO utils.Index, dRhoF, dRhoUF, dEnerF *utils.Matrix) {
	// Periodic Boundary condition
	fRho, fRhoU, fEner := c.RoeFlux(Rho, RhoU, Ener, RhoF, RhoUF, EnerF, vmapI, vmapO)
	dRhoF.AssignVector(c.El.MapI, fRho)
	dRhoUF.AssignVector(c.El.MapI, fRhoU)
	dEnerF.AssignVector(c.El.MapI, fEner)
	dRhoF.AssignVector(c.El.MapO, fRho)
	dRhoUF.AssignVector(c.El.MapO, fRhoU)
	dEnerF.AssignVector(c.El.MapO, fEner)
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

func (c *Euler) AveFlux(Rho, RhoU, Ener, RhoF, RhoUF, EnerF utils.Matrix, vmapM, vmapP utils.Index) (fRho, fRhoU, fEner utils.Matrix) {
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

func (c *Euler) LaxFlux(Rho, RhoU, Ener, RhoF, RhoUF, EnerF utils.Matrix, vmapM, vmapP utils.Index) (fRho, fRhoU, fEner utils.Matrix) {
	var (
		el       = c.El
		s        = c.State
		nrF, ncF = el.Nfp * el.NFaces, el.K
	)
	// Compute Lax-Friedrichs flux
	// Face jumps in primary and flux variables
	fJump := func(U utils.Matrix) (dU utils.Matrix) {
		dU = U.Subset(vmapM, nrF, ncF).Subtract(U.Subset(vmapP, nrF, ncF)).ElMul(el.NX)
		return
	}
	// Face flux average
	fAve := func(U utils.Matrix) (Uavg utils.Matrix) {
		Uavg = U.Subset(vmapM, nrF, ncF).Add(U.Subset(vmapP, nrF, ncF)).Scale(0.5)
		return
	}
	// Max eigenvalue
	LFc := s.LM.Subset(vmapM, nrF, ncF).Apply2(s.LM.Subset(vmapP, nrF, ncF), math.Max)
	// Compute numerical flux at faces
	fRho = fAve(RhoF).Add(fJump(Rho).ElMul(LFc).Scale(0.5))
	fRhoU = fAve(RhoUF).Add(fJump(RhoU).ElMul(LFc).Scale(0.5))
	fEner = fAve(EnerF).Add(fJump(Ener).ElMul(LFc).Scale(0.5))
	return
}

func (c *Euler) RoeFlux(Rho, RhoU, Ener, RhoF, RhoUF, EnerF utils.Matrix, vmapM, vmapP utils.Index) (fRho, fRhoU, fEner utils.Matrix) {
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
		Ul = U.Subset(vmapM, nrF, ncF)
		return
	}
	fR := func(U utils.Matrix) (Ul utils.Matrix) {
		Ul = U.Subset(vmapP, nrF, ncF)
		return
	}
	// Face jumps in primary and flux variables
	fJump := func(U utils.Matrix) (dU utils.Matrix) {
		dU = U.Subset(vmapP, nrF, ncF).Subtract(U.Subset(vmapM, nrF, ncF))
		return
	}
	// Face average
	fAve := func(U utils.Matrix) (Uavg utils.Matrix) {
		Uavg = U.Subset(vmapP, nrF, ncF).Add(U.Subset(vmapM, nrF, ncF)).Scale(0.5)
		return
	}
	//	Calculate the Roe Averaged variables
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

func (fs *FieldState) Update(Rho, RhoU, Ener utils.Matrix, c *Euler) (RhoFull, RhoUFull, EnerFull, RhoF, RhoUF, EnerF utils.Matrix) {
	//	Calorically Perfect Gas, R = 1
	var (
		Gamma                  = 1.4
		Cv                     = 1. / (Gamma - 1.)
		Cp                     = Gamma * Cv
		FluxRanger, FluxSubset = c.FluxRanger, c.FluxSubset
		limiter                = c.useLimiter
		slopeLimiterM          = 20.
		el                     = c.El
	)
	// Copy the solution at solution points, then interpolate to the edges of the Flux space
	RhoFull = c.InterpolateBoundaries(utils.NewMatrix(FluxRanger.Dims()).AssignVector(FluxSubset, Rho))
	RhoUFull = c.InterpolateBoundaries(utils.NewMatrix(FluxRanger.Dims()).AssignVector(FluxSubset, RhoU))
	EnerFull = c.InterpolateBoundaries(utils.NewMatrix(FluxRanger.Dims()).AssignVector(FluxSubset, Ener))
	if limiter {
		// Slope Limit the solution fields
		RhoFull = el.SlopeLimitN(RhoFull, slopeLimiterM)
		RhoUFull = el.SlopeLimitN(RhoUFull, slopeLimiterM)
		EnerFull = el.SlopeLimitN(EnerFull, slopeLimiterM)
	}

	fs.U = RhoUFull.Copy().ElDiv(RhoFull) // Velocity
	fs.Q = fs.U.Copy().Apply2(RhoFull, func(u, rho float64) (q float64) {
		q = 0.5 * u * u * rho
		return
	})
	// (gamma-1.0)*(Ener - 0.5*(rhou).^2./rho)
	fs.Pres = EnerFull.Copy().Apply2(fs.Q, func(e, q float64) (p float64) {
		p = (e - q) * (fs.Gamma - 1)
		return
	})
	// sqrt(gamma*pres./rho)
	fs.CVel = fs.Pres.Copy().Apply2(RhoFull, func(p, rho float64) (cvel float64) {
		cvel = math.Sqrt(fs.Gamma * p / rho)
		return
	})
	//  abs(rhou./rho)+cvel
	fs.LM = fs.U.Copy().Apply2(fs.CVel, func(u, cvel float64) (lm float64) {
		lm = math.Abs(u) + cvel
		return
	})
	//	Temp = (Ener - 0.5*(rhou).^2./rho)./rho
	fs.Temp = EnerFull.Copy().Apply3(fs.Q, RhoFull, func(e, q, rho float64) (temp float64) {
		temp = (e - q) / rho
		return
	})
	fs.Ht = fs.Q.Copy().Apply3(RhoFull, fs.Temp, func(q, rho, temp float64) (ht float64) {
		ht = Cp * (q/rho + temp)
		return
	})
	RhoF = RhoUFull
	RhoUF = fs.Q.Copy().Apply2(fs.Pres, func(q, pres float64) (rhouf float64) {
		rhouf = 2*q + pres
		return
	})
	EnerF = EnerFull.Copy().Apply3(fs.Pres, fs.U, func(ener, pres, u float64) (enerf float64) {
		enerf = u * (ener + pres)
		return
	})
	fs.U.SetReadOnly("U")
	fs.Q.SetReadOnly("Q")
	fs.Pres.SetReadOnly("Pres")
	fs.CVel.SetReadOnly("CVel")
	fs.LM.SetReadOnly("LM")
	fs.Ht.SetReadOnly("Ht")
	fs.Temp.SetReadOnly("Temp")
	return
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
*/

package Euler2D

import (
	"fmt"
	"math"
	"runtime"
	"time"

	"github.com/notargets/gocfd/types"

	"github.com/notargets/gocfd/DG2D"

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
	MeshFile          string
	CFL, FinalTime    float64
	FS                *FreeStream
	dfr               *DG2D.DFR2D
	chart             ChartState
	FluxCalcAlgo      FluxType
	Case              InitType
	AnalyticSolution  ExactState
	FluxCalcMock      func(Q [4]float64) (Fx, Fy [4]float64) // For testing
	SortedEdgeKeys    EdgeKeySlice
	ParallelDegree    int // Number of go routines to use for parallel execution
	Partitions        *PartitionMap
	LocalTimeStepping bool
	MaxIterations     int
	// Below are partitioned by K (elements) in the first slice
	Q                    [][4]utils.Matrix // Solution variables, stored at solution point locations, Np_solution x K
	SolutionX, SolutionY []utils.Matrix
}

func NewEuler(FinalTime float64, N int, meshFile string, CFL float64,
	fluxType FluxType, Case InitType, ProcLimit int,
	Minf, Gamma, Alpha float64, LocalTime bool, MaxIterations int, plotMesh, verbose bool) (c *Euler) {
	c = &Euler{
		MeshFile:          meshFile,
		CFL:               CFL,
		FinalTime:         FinalTime,
		FluxCalcAlgo:      fluxType,
		Case:              Case,
		LocalTimeStepping: LocalTime,
		MaxIterations:     MaxIterations,
		FS:                NewFreeStream(Minf, Gamma, Alpha),
	}
	c.FluxCalcMock = c.FluxCalc
	if ProcLimit != 0 {
		c.ParallelDegree = ProcLimit
	} else {
		c.ParallelDegree = runtime.NumCPU()
	}
	runtime.GOMAXPROCS(runtime.NumCPU())

	if len(meshFile) == 0 {
		return
	}

	c.dfr = DG2D.NewDFR2D(N, plotMesh, meshFile)
	if c.ParallelDegree > c.dfr.K {
		c.ParallelDegree = 1
	}
	c.Partitions = NewPartitionMap(c.ParallelDegree, c.dfr.K)

	if verbose {
		fmt.Printf("Euler Equations in 2 Dimensions\n")
		fmt.Printf("Using %d go routines in parallel\n", c.Partitions.ParallelDegree)
		fmt.Printf("Solving %s\n", c.Case.Print())
		if c.Case == FREESTREAM {
			fmt.Printf("Mach Infinity = %8.5f, Angle of Attack = %8.5f\n", Minf, Alpha)
		}
	}

	// Setup the key for edge calculations, useful for parallelizing the process
	c.SortedEdgeKeys = make(EdgeKeySlice, len(c.dfr.Tris.Edges))
	var i int
	for en := range c.dfr.Tris.Edges {
		c.SortedEdgeKeys[i] = en
		i++
	}
	c.SortedEdgeKeys.Sort()

	// Initialize solution
	c.Q = make([][4]utils.Matrix, c.Partitions.ParallelDegree) // Allocate shards to store solution
	switch c.Case {
	case FREESTREAM:
		NP := c.Partitions.ParallelDegree
		for np := 0; np < NP; np++ {
			Kmax := c.Partitions.GetBucketDimension(np)
			c.Q[np] = c.InitializeFS(Kmax)
		}
	case IVORTEX:
		c.FS.Qinf = [4]float64{1, 1, 0, 3}
		c.FS.Pinf = c.FS.GetFlowFunction(c.FS.Qinf, StaticPressure)
		c.FS.QQinf = c.FS.GetFlowFunction(c.FS.Qinf, DynamicPressure)
		c.FS.Cinf = c.FS.GetFlowFunction(c.FS.Qinf, SoundSpeed)
		c.SolutionX = c.ShardByK(c.dfr.SolutionX)
		c.SolutionY = c.ShardByK(c.dfr.SolutionY)
		NP := c.Partitions.ParallelDegree
		for np := 0; np < NP; np++ {
			c.AnalyticSolution, c.Q[np] = c.InitializeIVortex(c.SolutionX[np], c.SolutionY[np])
		}
		// Set "Wall" BCs to IVortex
		var count int
		for _, e := range c.dfr.Tris.Edges {
			if e.BCType == types.BC_Wall {
				count++
				e.BCType = types.BC_IVortex
			}
		}
		if verbose {
			fmt.Printf("\tReplaced %d Wall boundary conditions with analytic BC_IVortex\n", count)
		}
	default:
		panic("unknown case type")
	}
	if verbose {
		fmt.Printf("Algorithm: %s\n", c.FluxCalcAlgo.Print())
		fmt.Printf("CFL = %8.4f, Polynomial Degree N = %d (1 is linear), Num Elements K = %d\n\n\n", CFL, N, c.dfr.K)
	}
	return
}

func (c *Euler) Solve(pm *PlotMeta) {
	var (
		FinalTime            = c.FinalTime
		Time, dt             float64
		steps                int
		finished             bool
		Np                   = c.dfr.SolutionElement.Np
		plotQ                = pm.Plot
		pts                  = c.Partitions
		NPar                 = pts.ParallelDegree
		Q1, Q2, Q3, Residual [][4]utils.Matrix
		F_RT_DOF             [][4]utils.Matrix
		DT                   []utils.Matrix
	)
	c.PrintInitialization(FinalTime)
	// Initialize memory for RHS
	Q1, Q2, Q3, Residual, F_RT_DOF = make([][4]utils.Matrix, NPar), make([][4]utils.Matrix, NPar),
		make([][4]utils.Matrix, NPar), make([][4]utils.Matrix, NPar), make([][4]utils.Matrix, NPar)
	DT = make([]utils.Matrix, NPar)
	Jdet, Jinv := c.ShardByK(c.dfr.Jdet), c.ShardByK(c.dfr.Jinv)

	for np := 0; np < NPar; np++ {
		Kmax := pts.GetBucketDimension(np)
		for n := 0; n < 4; n++ {
			Q1[np][n] = utils.NewMatrix(Np, Kmax)
			Q2[np][n] = utils.NewMatrix(Np, Kmax)
			Q3[np][n] = utils.NewMatrix(Np, Kmax)
			F_RT_DOF[np][n] = utils.NewMatrix(Np, Kmax)
		}
		DT[np] = utils.NewMatrix(Np, Kmax)
	}
	start := time.Now()
	for !finished {
		Residual, dt = c.RungeKutta4SSP(Time, Jdet, Jinv, DT, F_RT_DOF, c.Q, Q1, Q2, Q3)
		steps++
		Time += dt
		finished = c.CheckIfFinished(Time, FinalTime, steps)
		if finished || steps%pm.StepsBeforePlot == 0 || steps == 1 {
			QQ := c.RecombineShardsKBy4(c.Q)
			ResidFull := c.RecombineShardsKBy4(Residual)
			c.PrintUpdate(Time, dt, steps, QQ, ResidFull, plotQ, pm)
		}
	}
	elapsed := time.Now().Sub(start)
	c.PrintFinal(elapsed, steps)
}

func (c *Euler) RungeKutta4SSP(Time float64, Jdet, Jinv []utils.Matrix,
	DT []utils.Matrix, F_RT_DOF, Q0, Q1, Q2, Q3 [][4]utils.Matrix) (Residual [][4]utils.Matrix, dt float64) {
	var (
		Np     = c.dfr.SolutionElement.Np
		dT     float64
		NP     = c.Partitions.ParallelDegree
		Q_Face [][4]utils.Matrix
	)
	Q_Face, Residual = make([][4]utils.Matrix, NP), make([][4]utils.Matrix, NP)
	for np := 0; np < NP; np++ {
		for n := 0; n < 4; n++ {
			Residual[np][n] = Q1[np][n] // optimize memory using an alias
		}
		Kmax := c.Partitions.GetBucketDimension(np)
		Q_Face[np] = c.PrepareEdgeFlux(Kmax, Jdet[np], Jinv[np], F_RT_DOF[np], Q0[np], Time)
	}
	dt = c.CalculateDT(DT, Q_Face) // Operates on edges/ internally parallel
	if Time+dt > c.FinalTime {
		dt = c.FinalTime - Time
	}
	c.ParallelSetNormalFluxOnEdges(Time, F_RT_DOF, Q_Face) // Must sync parallel before calling
	for np := 0; np < NP; np++ {
		Kmax := c.Partitions.GetBucketDimension(np)
		qD := Get4DP(Q0[np])
		q1D := Get4DP(Q1[np])
		dtD := DT[np].Data()
		rhsQ := c.RHS(Kmax, Jdet[np], Q_Face[np], F_RT_DOF[np], Q0[np], Time)
		rhsD := Get4DP(rhsQ)
		for n := 0; n < 4; n++ {
			for i := 0; i < Kmax*Np; i++ {
				if c.LocalTimeStepping {
					dT = dtD[i]
				}
				q1D[n][i] = qD[n][i] + 0.5*rhsD[n][i]*dT
			}
		}
		Q_Face[np] = c.PrepareEdgeFlux(Kmax, Jdet[np], Jinv[np], F_RT_DOF[np], Q1[np], Time)
	}
	c.ParallelSetNormalFluxOnEdges(Time, F_RT_DOF, Q_Face) // Must sync parallel before calling
	for np := 0; np < NP; np++ {
		Kmax := c.Partitions.GetBucketDimension(np)
		q1D := Get4DP(Q1[np])
		q2D := Get4DP(Q2[np])
		dtD := DT[np].Data()
		rhsQ := c.RHS(Kmax, Jdet[np], Q_Face[np], F_RT_DOF[np], Q1[np], Time)
		rhsD := Get4DP(rhsQ)
		for n := 0; n < 4; n++ {
			for i := 0; i < Kmax*Np; i++ {
				if c.LocalTimeStepping {
					dT = dtD[i]
				}
				q2D[n][i] = q1D[n][i] + 0.25*rhsD[n][i]*dT
			}
		}
		Q_Face[np] = c.PrepareEdgeFlux(Kmax, Jdet[np], Jinv[np], F_RT_DOF[np], Q2[np], Time)
	}
	c.ParallelSetNormalFluxOnEdges(Time, F_RT_DOF, Q_Face) // Must sync parallel before calling
	for np := 0; np < NP; np++ {
		Kmax := c.Partitions.GetBucketDimension(np)
		qD := Get4DP(Q0[np])
		q2D := Get4DP(Q2[np])
		q3D := Get4DP(Q3[np])
		dtD := DT[np].Data()
		rhsQ := c.RHS(Kmax, Jdet[np], Q_Face[np], F_RT_DOF[np], Q2[np], Time)
		rhsD := Get4DP(rhsQ)
		for n := 0; n < 4; n++ {
			for i := 0; i < Kmax*Np; i++ {
				if c.LocalTimeStepping {
					dT = dtD[i]
				}
				q3D[n][i] = (1. / 3.) * (2*qD[n][i] + q2D[n][i] + rhsD[n][i]*dT)
			}
		}
		Q_Face[np] = c.PrepareEdgeFlux(Kmax, Jdet[np], Jinv[np], F_RT_DOF[np], Q3[np], Time)
	}
	c.ParallelSetNormalFluxOnEdges(Time, F_RT_DOF, Q_Face) // Must sync parallel before calling
	for np := 0; np < NP; np++ {
		Kmax := c.Partitions.GetBucketDimension(np)
		qD := Get4DP(Q0[np])
		q3D := Get4DP(Q3[np])
		resD := Get4DP(Residual[np])
		dtD := DT[np].Data()
		rhsQ := c.RHS(Kmax, Jdet[np], Q_Face[np], F_RT_DOF[np], Q3[np], Time)
		rhsD := Get4DP(rhsQ)
		for n := 0; n < 4; n++ {
			for i := 0; i < Kmax*Np; i++ {
				if c.LocalTimeStepping {
					dT = dtD[i]
				}
				resD[n][i] = q3D[n][i] + 0.25*rhsD[n][i]*dT - qD[n][i]
				qD[n][i] += resD[n][i]
			}
		}
	}
	return
}

func (c *Euler) RHS(Kmax int, Jdet utils.Matrix, Q_Face, F_RT_DOF, Q [4]utils.Matrix, Time float64) (RHSCalc [4]utils.Matrix) {
	/*
				Calculate the RHS of the equation:
				dQ/dt = -div(F,G)
				Where:
					Q = [rho, rhoU, rhoV, E]
					F = [rhoU, rhoU*u+p, rhoV*u, u*(E+p)]
					G = [rhoV, rhoU*v, rhoV*v+p, v*(E+p)]

		    	The divergence div(F,G) is calculated using a Raviart Thomas finite element with flux (F,G) values on the faces
				of the element "injected" via calculation of a physical flux on those faces, and the (F,G) values in the interior
				of the element taken directly from the solution values (Q).
	*/
	for n := 0; n < 4; n++ {
		RHSCalc[n] = c.dfr.FluxElement.DivInt.Mul(F_RT_DOF[n]) // Calculate divergence for the internal node points
		c.DivideByJacobian(Kmax, c.dfr.FluxElement.Nint, Jdet, RHSCalc[n].Data(), -1)
	}
	return
}

func (c *Euler) DivideByJacobian(Kmax, Imax int, Jdet utils.Matrix, data []float64, scale float64) {
	var (
		JdetD = Jdet.Data()
	)
	for k := 0; k < Kmax; k++ {
		for i := 0; i < Imax; i++ {
			ind := k + i*Kmax
			data[ind] /= JdetD[k] * scale // Multiply divergence by -1 to produce the RHS
		}
	}
}

func (c *Euler) PrepareEdgeFlux(Kmax int, Jdet, Jinv utils.Matrix, F_RT_DOF, Q [4]utils.Matrix, Time float64) (Q_Face [4]utils.Matrix) {
	var (
		Np    = c.dfr.FluxElement.Np
		fdofD = Get4DP(F_RT_DOF)
	)
	/*
		Solver approach:
		0) Solution is stored on sol points as Q
		0a) Flux is computed and stored in X, Y component projections in the 2*Nint front of F_RT_DOF
		1) Solution is extrapolated to edge points in Q_Face from Q
		2) Edges are traversed, flux is calculated and projected onto edge face normals, scaled and placed into F_RT_DOF
	*/
	/*
		Zero out DOF storage to promote easier bug avoidance
	*/
	for n := 0; n < 4; n++ {
		for i := 0; i < Kmax*Np; i++ {
			fdofD[n][i] = 0.
		}
	}
	c.SetNormalFluxInternal(Kmax, Jdet, Jinv, F_RT_DOF, Q) // Updates F_RT_DOF with values from Q
	Q_Face = c.InterpolateSolutionToEdges(Q)               // Interpolates Q_Face values from Q
	return
}

func (c *Euler) SetNormalFluxInternal(Kmax int, Jdet, Jinv utils.Matrix, F_RT_DOF, Q [4]utils.Matrix) {
	var (
		Nint  = c.dfr.FluxElement.Nint
		qD    = Get4DP(Q)
		fdofD = Get4DP(F_RT_DOF)
	)
	// Calculate flux and project into R and S (transformed) directions for the internal points
	for k := 0; k < Kmax; k++ {
		for i := 0; i < Nint; i++ {
			ind := k + i*Kmax
			ind2 := k + (i+Nint)*Kmax
			Fr, Fs := c.CalculateFluxTransformed(k, Kmax, i, Jdet, Jinv, qD)
			for n := 0; n < 4; n++ {
				fdofD[n][ind], fdofD[n][ind2] = Fr[n], Fs[n]
			}
		}
	}
}

func (c *Euler) CheckIfFinished(Time, FinalTime float64, steps int) (finished bool) {
	if Time >= FinalTime || steps >= c.MaxIterations {
		finished = true
	}
	return
}
func (c *Euler) PrintInitialization(FinalTime float64) {
	fmt.Printf("Using mesh from file: [%s]\n", c.MeshFile)
	if c.LocalTimeStepping {
		fmt.Printf("Solving until Max Iterations = %d\n", c.MaxIterations)
		fmt.Printf("    iter                ")
	} else {
		fmt.Printf("Solving until finaltime = %8.5f\n", FinalTime)
		fmt.Printf("    iter    time  min_dt")
	}
	fmt.Printf("       Res0       Res1       Res2")
	fmt.Printf("       Res3         L1         L2\n")
}
func (c *Euler) PrintUpdate(Time, dt float64, steps int, Q, Residual [4]utils.Matrix, plotQ bool, pm *PlotMeta) {
	format := "%11.4e"
	if plotQ {
		c.PlotQ(Q, pm) // wait till we implement time iterative frame updates
	}
	if c.LocalTimeStepping {
		fmt.Printf("%10d              ", steps)
	} else {
		fmt.Printf("%8d%8.5f%8.5f", steps, Time, dt)
	}
	var l1, l2 float64
	for n := 0; n < 4; n++ {
		maxR := Residual[n].Max()
		fmt.Printf(format, maxR)
		if maxR > l1 {
			l1 = maxR
		}
		l2 += maxR * maxR
	}
	fmt.Printf(format, l1)
	fmt.Printf(format, math.Sqrt(l2)/4.)
	fmt.Printf("\n")
}
func (c *Euler) PrintFinal(elapsed time.Duration, steps int) {
	rate := float64(elapsed.Microseconds()) / (float64(c.dfr.K * steps))
	fmt.Printf("\nRate of execution = %8.5f us/(element*iteration) over %d iterations\n", rate, steps)
}

func (c *Euler) ShardByK(A utils.Matrix) (pA []utils.Matrix) {
	var (
		NP      = c.Partitions.ParallelDegree
		_, Imax = A.Dims()
		aD      = A.Data()
	)
	pA = make([]utils.Matrix, NP)
	for np := 0; np < NP; np++ {
		kMin, kMax := c.Partitions.GetBucketRange(np)
		Kmax := c.Partitions.GetBucketDimension(np)
		pA[np] = utils.NewMatrix(Kmax, Imax)
		pAD := pA[np].Data()
		for k := kMin; k < kMax; k++ {
			pk := k - kMin
			for i := 0; i < Imax; i++ {
				ind := k + i*c.dfr.K
				pind := pk + Kmax*i
				pAD[pind] = aD[ind]
			}
		}
	}
	return
}

func (c *Euler) RecombineShardsK(pA []utils.Matrix) (A utils.Matrix) {
	var (
		NP      = c.Partitions.ParallelDegree
		_, Imax = pA[0].Dims()
	)
	A = utils.NewMatrix(c.dfr.K, Imax)
	aD := A.Data()
	for np := 0; np < NP; np++ {
		kMin, kMax := c.Partitions.GetBucketRange(np)
		Kmax := c.Partitions.GetBucketDimension(np)
		pAD := pA[np].Data()
		for k := kMin; k < kMax; k++ {
			pk := k - kMin
			for i := 0; i < Imax; i++ {
				ind := k + i*c.dfr.K
				pind := pk + Kmax*i
				aD[ind] = pAD[pind]
			}
		}
	}
	return
}

func (c *Euler) RecombineShardsKBy4(pA [][4]utils.Matrix) (A [4]utils.Matrix) {
	var (
		NP      = c.Partitions.ParallelDegree
		_, Imax = pA[0][0].Dims()
	)
	for np := 0; np < NP; np++ {
		for n := 0; n < 4; n++ {
			A[n] = utils.NewMatrix(c.dfr.K, Imax)
			aD := A[n].Data()
			kMin, kMax := c.Partitions.GetBucketRange(np)
			Kmax := c.Partitions.GetBucketDimension(np)
			pAD := pA[np][n].Data()
			for k := kMin; k < kMax; k++ {
				pk := k - kMin
				for i := 0; i < Imax; i++ {
					ind := k + i*c.dfr.K
					pind := pk + Kmax*i
					aD[ind] = pAD[pind]
				}
			}
		}
	}
	return
}

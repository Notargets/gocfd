package Euler2D

import (
	"fmt"
	"math"
	"syscall"
	"time"

	"github.com/notargets/gocfd/types"

	"github.com/notargets/gocfd/DG2D"

	"github.com/notargets/gocfd/utils"

	"github.com/pkg/profile"
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
	profile           bool // Generate a CPU profile of the solver
	FluxCalcAlgo      FluxType
	Case              InitType
	AnalyticSolution  ExactState
	FluxCalcMock      func(rho, rhoU, rhoV, E float64) (Fx, Fy [4]float64) // For testing
	SortedEdgeKeys    []EdgeKeySlice                                       // Buckets, one for each parallel partition
	Partitions        *PartitionMap                                        // mapping of elements into bins for parallelism
	LocalTimeStepping bool
	MaxIterations     int
	// Below are partitioned by K (elements) in the first slice
	Q                    [][4]utils.Matrix // Solution variables, stored at solution point locations, Np_solution x K
	SolutionX, SolutionY []utils.Matrix
}

func NewEuler(FinalTime float64, N int, meshFile string, CFL float64, fluxType FluxType, Case InitType,
	ProcLimit int, Minf, Gamma, Alpha float64, LocalTime bool,
	MaxIterations int, plotMesh, verbose, profile bool) (c *Euler) {
	c = &Euler{
		MeshFile:          meshFile,
		CFL:               CFL,
		FinalTime:         FinalTime,
		FluxCalcAlgo:      fluxType,
		Case:              Case,
		LocalTimeStepping: LocalTime,
		MaxIterations:     MaxIterations,
		FS:                NewFreeStream(Minf, Gamma, Alpha),
		profile:           profile,
	}
	c.FluxCalcMock = c.FluxCalcBase

	if len(meshFile) == 0 {
		return
	}

	// Read mesh file, initialize geometry and finite elements
	c.dfr = DG2D.NewDFR2D(N, plotMesh, meshFile)

	c.SetParallelDegree(ProcLimit, c.dfr.K) // Must occur after determining the number of elements

	c.PartitionEdgesByK() // Setup the key for edge calculations, useful for parallelizing the process

	c.InitializeSolution(verbose)

	if verbose {
		fmt.Printf("Euler Equations in 2 Dimensions\n")
		fmt.Printf("Using %d go routines in parallel\n", c.Partitions.ParallelDegree)
		fmt.Printf("Solving %s\n", c.Case.Print())
		if c.Case == FREESTREAM {
			fmt.Printf("Mach Infinity = %8.5f, Angle of Attack = %8.5f\n", Minf, Alpha)
		}
		fmt.Printf("Algorithm: %s\n", c.FluxCalcAlgo.Print())
		fmt.Printf("CFL = %8.4f, Polynomial Degree N = %d (1 is linear), Num Elements K = %d\n\n\n", CFL, N, c.dfr.K)
	}
	return
}

func (c *Euler) Solve(pm *PlotMeta) {
	var (
		FinalTime = c.FinalTime
		steps     int
		finished  bool
		plotQ     = pm.Plot
	)
	if c.profile {
		defer profile.Start().Stop()
	}

	c.PrintInitialization(FinalTime)

	rk := c.NewRungeKuttaSSP()

	elapsed := time.Duration(0)
	rk.StartWorkers(c)
	var start time.Time
	for !finished {
		start = time.Now()
		rk.Step(c)
		elapsed += time.Now().Sub(start)
		steps++
		rk.Time += rk.GlobalDT
		finished = c.CheckIfFinished(rk.Time, FinalTime, steps)
		if finished || steps%pm.StepsBeforePlot == 0 || steps == 1 {
			c.PrintUpdate(rk.Time, rk.GlobalDT, steps, c.Q, rk.Residual, plotQ, pm)
		}
	}
	c.PrintFinal(elapsed, steps)
}

type RungeKutta4SSP struct {
	Jdet, Jinv        []utils.Matrix    // Sharded mesh Jacobian and inverse transform
	RHSQ, Q_Face      [][4]utils.Matrix // State used for matrix multiplies within the time step algorithm
	Q1, Q2, Q3        [][4]utils.Matrix // Intermediate solution state
	Residual          [][4]utils.Matrix // Used for reporting, aliased to Q1
	F_RT_DOF          [][4]utils.Matrix // Normal flux used for divergence
	DT                []utils.Matrix    // Local time step storage
	MaxWaveSpeed      []float64         // Shard max wavespeed
	GlobalDT, Time    float64
	Kmax              []int // Local element count (dimension: Kmax[ParallelDegree])
	Np, Nedge, NpFlux int   // Number of points in solution, edge and flux total
	toWorker          []chan struct{}
	fromWorkers       chan int8
}

func (c *Euler) NewRungeKuttaSSP() (rk *RungeKutta4SSP) {
	var (
		pm   = c.Partitions
		NPar = pm.ParallelDegree
	)
	rk = &RungeKutta4SSP{
		Jdet:         c.ShardByKTranspose(c.dfr.Jdet),
		Jinv:         c.ShardByKTranspose(c.dfr.Jinv),
		RHSQ:         make([][4]utils.Matrix, NPar),
		Q_Face:       make([][4]utils.Matrix, NPar),
		Q1:           make([][4]utils.Matrix, NPar),
		Q2:           make([][4]utils.Matrix, NPar),
		Q3:           make([][4]utils.Matrix, NPar),
		Residual:     make([][4]utils.Matrix, NPar),
		F_RT_DOF:     make([][4]utils.Matrix, NPar),
		DT:           make([]utils.Matrix, NPar),
		MaxWaveSpeed: make([]float64, NPar),
		Kmax:         make([]int, NPar),
		Np:           c.dfr.SolutionElement.Np,
		Nedge:        c.dfr.FluxElement.Nedge,
		NpFlux:       c.dfr.FluxElement.Np,
		fromWorkers:  make(chan int8, NPar),
		toWorker:     make([]chan struct{}, NPar),
	}
	for np := 0; np < NPar; np++ {
		rk.toWorker[np] = make(chan struct{}, 1)
	}
	// Initialize memory for RHS
	for np := 0; np < NPar; np++ {
		rk.Kmax[np] = pm.GetBucketDimension(np)
		for n := 0; n < 4; n++ {
			rk.Q1[np][n] = utils.NewMatrix(rk.Np, rk.Kmax[np])
			rk.Q2[np][n] = utils.NewMatrix(rk.Np, rk.Kmax[np])
			rk.Q3[np][n] = utils.NewMatrix(rk.Np, rk.Kmax[np])
			rk.RHSQ[np][n] = utils.NewMatrix(rk.Np, rk.Kmax[np])
			rk.F_RT_DOF[np][n] = utils.NewMatrix(rk.NpFlux, rk.Kmax[np])
			rk.Q_Face[np][n] = utils.NewMatrix(rk.Nedge*3, rk.Kmax[np])
		}
		rk.DT[np] = utils.NewMatrix(rk.Np, rk.Kmax[np])
	}
	rk.Residual = rk.Q1
	return
}

func (rk *RungeKutta4SSP) StartWorkers(c *Euler) {
	var (
		pm = c.Partitions
		NP = pm.ParallelDegree
	)
	for np := 0; np < NP; np++ {
		go rk.StepWorker(c, np, rk.toWorker[np], rk.fromWorkers)
	}
}

func (rk *RungeKutta4SSP) Step(c *Euler) {
	/*
		This is the controller thread - it manages and synchronizes the workers
		NOTE: Any CPU work done in this thread makes the workers wait - be careful!
	*/
	var (
		pm = c.Partitions
		NP = pm.ParallelDegree
	)
	for {
		// Advance the workers one step by sending the "go" message. They are blocked until they receive this.
		for np := 0; np < NP; np++ {
			rk.toWorker[np] <- struct{}{}
		}
		currentStep := int8(-1)
		var ts int8
		for np := 0; np < NP; np++ {
			// Each worker sends it's completed step number, we check to make sure they are all the same
			ts = <-rk.fromWorkers
			if currentStep == -1 {
				currentStep = ts
			}
			if ts != currentStep {
				err := fmt.Errorf("[%d]incorrect state, ts = %d, currentStep = %d\n", np, ts, currentStep)
				panic(err)
			}
		}
		// Workers are blocked below here - make sure no significant CPU work is done here unless abs necessary!
		switch {
		case currentStep == 2:
			// After step 2, we need to consolidate the local max wavespeeds to calculate global dt
			if !c.LocalTimeStepping {
				rk.calculateGlobalDT(c)
			}
		case currentStep < 0:
			return
		}
	}
}

func (rk *RungeKutta4SSP) WorkerDone(subStepP *int8, fromWorker chan int8, isDone bool) {
	if isDone {
		*subStepP = -1
	} else {
		*subStepP++
	}
	fromWorker <- *subStepP // Inform controller what step we've just completed
}

func (rk *RungeKutta4SSP) StepWorker(c *Euler, myThread int, fromController chan struct{}, toController chan int8) {
	var (
		Q0                           = c.Q
		Np                           = rk.Np
		Kmax, Jdet, Jinv, F_RT_DOF   = rk.Kmax, rk.Jdet, rk.Jinv, rk.F_RT_DOF
		DT, Q_Face, Q1, Q2, Q3, RHSQ = rk.DT, rk.Q_Face, rk.Q1, rk.Q2, rk.Q3, rk.RHSQ
		dT                           float64
		subStep                      int8
	)
	for {
		_ = <-fromController // Block until parent sends "go"
		if c.LocalTimeStepping {
			// Setup local time stepping
			for k := 0; k < Kmax[myThread]; k++ {
				DT[myThread].DataP[k] = -100 // Global
			}
		}
		c.PrepareEdgeFlux(Kmax[myThread], Jdet[myThread], Jinv[myThread], F_RT_DOF[myThread], Q0[myThread], Q_Face[myThread])
		rk.WorkerDone(&subStep, toController, false)

		_ = <-fromController // Block until parent sends "go"
		if !c.LocalTimeStepping {
			rk.MaxWaveSpeed[myThread] = c.SetNormalFluxOnEdges(rk.Time, true, Jdet, DT, F_RT_DOF, Q_Face, c.SortedEdgeKeys[myThread]) // Global
		} else {
			c.SetNormalFluxOnEdges(rk.Time, true, Jdet, DT, F_RT_DOF, Q_Face, c.SortedEdgeKeys[myThread]) // Global
		}
		rk.WorkerDone(&subStep, toController, false)

		_ = <-fromController // Block until parent sends "go"
		if c.LocalTimeStepping {
			// Replicate local time step to the other solution points for each k
			for k := 0; k < Kmax[myThread]; k++ {
				DT[myThread].DataP[k] = c.CFL / DT[myThread].DataP[k]
			}
			for i := 1; i < c.dfr.SolutionElement.Np; i++ {
				for k := 0; k < Kmax[myThread]; k++ {
					ind := k + Kmax[myThread]*i
					DT[myThread].DataP[ind] = DT[myThread].DataP[k]
				}
			}
		}
		c.RHS(Kmax[myThread], Jdet[myThread], F_RT_DOF[myThread], RHSQ[myThread])
		dT = rk.GlobalDT
		for n := 0; n < 4; n++ {
			for i := 0; i < Kmax[myThread]*Np; i++ {
				if c.LocalTimeStepping {
					dT = DT[myThread].DataP[i]
				}
				Q1[myThread][n].DataP[i] = Q0[myThread][n].DataP[i] + 0.5*RHSQ[myThread][n].DataP[i]*dT
			}
		}
		c.PrepareEdgeFlux(Kmax[myThread], Jdet[myThread], Jinv[myThread], F_RT_DOF[myThread], Q1[myThread], Q_Face[myThread])
		rk.WorkerDone(&subStep, toController, false)

		_ = <-fromController                                                                           // Block until parent sends "go"
		c.SetNormalFluxOnEdges(rk.Time, false, Jdet, DT, F_RT_DOF, Q_Face, c.SortedEdgeKeys[myThread]) // Global
		rk.WorkerDone(&subStep, toController, false)

		_ = <-fromController // Block until parent sends "go"
		c.RHS(Kmax[myThread], Jdet[myThread], F_RT_DOF[myThread], RHSQ[myThread])
		dT = rk.GlobalDT
		for n := 0; n < 4; n++ {
			for i := 0; i < Kmax[myThread]*Np; i++ {
				if c.LocalTimeStepping {
					dT = DT[myThread].DataP[i]
				}
				Q2[myThread][n].DataP[i] = Q1[myThread][n].DataP[i] + 0.25*RHSQ[myThread][n].DataP[i]*dT
			}
		}
		c.PrepareEdgeFlux(Kmax[myThread], Jdet[myThread], Jinv[myThread], F_RT_DOF[myThread], Q2[myThread], Q_Face[myThread])
		rk.WorkerDone(&subStep, toController, false)

		_ = <-fromController                                                                           // Block until parent sends "go"
		c.SetNormalFluxOnEdges(rk.Time, false, Jdet, DT, F_RT_DOF, Q_Face, c.SortedEdgeKeys[myThread]) // Global
		rk.WorkerDone(&subStep, toController, false)

		_ = <-fromController // Block until parent sends "go"
		c.RHS(Kmax[myThread], Jdet[myThread], F_RT_DOF[myThread], RHSQ[myThread])
		dT = rk.GlobalDT
		for n := 0; n < 4; n++ {
			for i := 0; i < Kmax[myThread]*Np; i++ {
				if c.LocalTimeStepping {
					dT = DT[myThread].DataP[i]
				}
				Q3[myThread][n].DataP[i] = (1. / 3.) * (2*Q0[myThread][n].DataP[i] + Q2[myThread][n].DataP[i] + RHSQ[myThread][n].DataP[i]*dT)
			}
		}
		c.PrepareEdgeFlux(Kmax[myThread], Jdet[myThread], Jinv[myThread], F_RT_DOF[myThread], Q3[myThread], Q_Face[myThread])
		rk.WorkerDone(&subStep, toController, false)

		_ = <-fromController                                                                           // Block until parent sends "go"
		c.SetNormalFluxOnEdges(rk.Time, false, Jdet, DT, F_RT_DOF, Q_Face, c.SortedEdgeKeys[myThread]) // Global
		rk.WorkerDone(&subStep, toController, false)

		_ = <-fromController // Block until parent sends "go"
		c.RHS(Kmax[myThread], Jdet[myThread], F_RT_DOF[myThread], RHSQ[myThread])
		// Note, we are re-using Q1 as storage for Residual here
		dT = rk.GlobalDT
		for n := 0; n < 4; n++ {
			for i := 0; i < Kmax[myThread]*Np; i++ {
				if c.LocalTimeStepping {
					dT = DT[myThread].DataP[i]
				}
				rk.Residual[myThread][n].DataP[i] = Q3[myThread][n].DataP[i] + 0.25*RHSQ[myThread][n].DataP[i]*dT - Q0[myThread][n].DataP[i]
				Q0[myThread][n].DataP[i] += rk.Residual[myThread][n].DataP[i]
			}
		}
		rk.WorkerDone(&subStep, toController, true)
	}
}

func (c *Euler) RHS(Kmax int, Jdet utils.Matrix, F_RT_DOF, RHSQ [4]utils.Matrix) {
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
		// Calculate divergence for the internal node points
		RHSQ[n] = c.dfr.FluxElement.DivInt.Mul(F_RT_DOF[n], RHSQ[n].DataP) // Re-use memory from RHSQ in the multiply
		c.DivideByJacobian(Kmax, c.dfr.FluxElement.Nint, Jdet, RHSQ[n].DataP, -1)
	}
	return
}

func (c *Euler) DivideByJacobian(Kmax, Imax int, Jdet utils.Matrix, data []float64, scale float64) {
	var (
		JdetD = Jdet.DataP
	)
	for k := 0; k < Kmax; k++ {
		for i := 0; i < Imax; i++ {
			ind := k + i*Kmax
			data[ind] /= JdetD[k] * scale // Multiply divergence by -1 to produce the RHS
		}
	}
}

func (c *Euler) PrepareEdgeFlux(Kmax int, Jdet, Jinv utils.Matrix, F_RT_DOF, Q, Q_Face [4]utils.Matrix) {
	var (
		Np    = c.dfr.FluxElement.Np
		fdofD = [4][]float64{F_RT_DOF[0].DataP, F_RT_DOF[1].DataP, F_RT_DOF[2].DataP, F_RT_DOF[3].DataP}
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
	c.InterpolateSolutionToEdges(Q, Q_Face)                // Interpolates Q_Face values from Q
	return
}

func (c *Euler) SetNormalFluxInternal(Kmax int, Jdet, Jinv utils.Matrix, F_RT_DOF, Q [4]utils.Matrix) {
	var (
		Nint  = c.dfr.FluxElement.Nint
		fdofD = [4][]float64{F_RT_DOF[0].DataP, F_RT_DOF[1].DataP, F_RT_DOF[2].DataP, F_RT_DOF[3].DataP}
	)
	// Calculate flux and project into R and S (transformed) directions for the internal points
	for k := 0; k < Kmax; k++ {
		for i := 0; i < Nint; i++ {
			ind := k + i*Kmax
			ind2 := k + (i+Nint)*Kmax
			Fr, Fs := c.CalculateFluxTransformed(k, Kmax, i, Jdet, Jinv, Q)
			for n := 0; n < 4; n++ {
				fdofD[n][ind], fdofD[n][ind2] = Fr[n], Fs[n]
			}
		}
	}
}

func (c *Euler) InitializeSolution(verbose bool) {
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
		c.FS.Pinf = c.FS.GetFlowFunctionQQ(c.FS.Qinf, StaticPressure)
		c.FS.QQinf = c.FS.GetFlowFunctionQQ(c.FS.Qinf, DynamicPressure)
		c.FS.Cinf = c.FS.GetFlowFunctionQQ(c.FS.Qinf, SoundSpeed)
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
func (c *Euler) PrintUpdate(Time, dt float64, steps int, Q, Residual [][4]utils.Matrix, plotQ bool, pm *PlotMeta) {
	format := "%11.4e"
	if plotQ {
		var QQ [4]utils.Matrix
		if plotQ {
			QQ = c.RecombineShardsKBy4(c.Q)
		} else {
			nilM := utils.Matrix{}
			QQ = [4]utils.Matrix{nilM, nilM, nilM, nilM}
		}
		c.PlotQ(QQ, pm) // wait till we implement time iterative frame updates
	}
	if c.LocalTimeStepping {
		fmt.Printf("%10d              ", steps)
	} else {
		fmt.Printf("%8d%8.5f%8.5f", steps, Time, dt)
	}
	var l1, l2 float64
	for n := 0; n < 4; n++ {
		var maxR float64
		for np := 0; np < c.Partitions.ParallelDegree; np++ {
			m := Residual[np][n].Max()
			if m > maxR {
				maxR = m
			}
		}
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
	instructionRate := 1344.282 // measured instructions per element per iteration
	iRate := instructionRate * float64(c.dfr.K*steps) / elapsed.Seconds()
	var mem syscall.Rusage
	_ = syscall.Getrusage(syscall.RUSAGE_SELF, &mem)
	usec, ssec := float64(mem.Utime.Nano())/1.e9, float64(mem.Stime.Nano())/1.e9
	fmt.Printf("OS Stats: Utime=%6.2f(s), Stime=%6.2f(s), RSS Memory=%6.2f(MB)\n", usec, ssec, float64(mem.Maxrss)/1024.)
	fmt.Printf("Estimated op rate = %8.5f Giga-ops / second, based on measure op rate of %8.5f per element per iteration\n",
		iRate/1000000000., instructionRate)
}

func (rk *RungeKutta4SSP) calculateGlobalDT(c *Euler) {
	var (
		pm = c.Partitions
		NP = pm.ParallelDegree
	)
	// Must be done in controller/sync process
	var wsMaxAll float64
	for np := 0; np < NP; np++ {
		wsMaxAll = math.Max(wsMaxAll, rk.MaxWaveSpeed[np])
	}
	rk.GlobalDT = c.CFL / wsMaxAll
	if rk.Time+rk.GlobalDT > c.FinalTime {
		rk.GlobalDT = c.FinalTime - rk.Time
	}
}

func (rk *RungeKutta4SSP) GetPreconditioner(c *Euler, rho, rhoU, rhoV, E float64) (P [4][4]float64) {
	/*
		[nx,ny] is the direction vector of the flux at the evaluated point:
			Flux = [F,G] = |[F,G]| * [nx,ny]
	*/
	var (
		gamma     = c.FS.Gamma
		sqrt, pow = math.Sqrt, utils.POW
		GM1       = gamma - 1.0
	)
	rhoU_2, rhoV_2 := rhoU*rhoU, rhoV*rhoV
	pressure := GM1 * (E - (rhoU_2+rhoV_2)/(rho*2.0))
	sig2 := -gamma*(rhoU_2) - gamma*(rhoV_2) + rhoU_2 + rhoV_2 + E*gamma*rho*2.0
	sig3 := 1.0 / rho
	sig3_2 := sig3 * sig3
	sig4 := 1.0 / sqrt(pow(pressure+rhoV_2*sig3, 2.0)+rhoU_2*rhoV_2*sig3_2)
	sig5 := 1.0 / sqrt(pow(pressure+rhoU_2*sig3, 2.0)+rhoU_2*rhoV_2*sig3_2)
	sig6 := rhoU_2 + rhoV_2
	sig7 := 1.0 / sqrt((sig2*sig2)*(sig3_2*sig3_2)*sig6)
	sig9 := pressure + rhoV_2*sig3
	sig10 := pressure + rhoU_2*sig3
	sig11 := rhoU_2 * rhoV_2 * (sig3 * sig3_2) * 2.0
	sig12 := sig2 + sig6*2.0 - gamma*sig6*2.0
	sig13 := rhoV_2 * sig3
	sig14 := rhoU_2 * sig3
	sig15 := 1.0 / sqrt(sig6)
	sig16 := ((sig3 * sig3) * sig6 * GM1) / 2.0
	sig17 := -GM1 * sig10
	sig18 := -GM1 * sig9
	sig19 := sig3 * sig3 * sig3 * sig3
	sig20 := sig2 * sig7 * sig12 * sig19
	P[0][1] = rhoU * sig15
	P[0][2] = rhoV * sig15
	P[1][0] = sig5 * (sig11 - sig10*(sig16-sig3*sig14)*2.0) * (-1.0 / 2.0)
	P[1][1] = rhoU * sig3 * sig5 * (sig10*2.0 + sig13 + sig17)
	P[1][2] = rhoV * sig3 * sig5 * (sig14 + sig17)
	P[1][3] = sig5 * GM1 * sig10
	P[2][0] = sig4 * (sig11 - sig9*(sig16-sig3*sig13)*2.0) * (-1.0 / 2.0)
	P[2][1] = rhoU * sig3 * sig4 * (sig13 + sig18)
	P[2][2] = rhoV * sig3 * sig4 * (sig9*2.0 + sig14 + sig18)
	P[2][3] = sig4 * GM1 * sig9
	P[3][0] = -sig2 * (sig3 * sig3 * sig3 * sig3 * sig3) * sig6 * sig7 * (sig2 - E*gamma*rho)
	P[3][1] = (rhoU * sig20) / 2.0
	P[3][2] = (rhoV * sig20) / 2.0
	P[3][3] = gamma * sig2 * (sig3 * sig3 * sig3) * sig6 * sig7
	return
}

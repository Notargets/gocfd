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
			var printMem bool
			if steps%100 == 0 {
				printMem = true
			}
			c.PrintUpdate(rk.Time, rk.GlobalDT, steps, c.Q, rk.Residual, plotQ, pm, printMem)
		}
	}
	c.PrintFinal(elapsed, steps)
}

type ElementImplicit struct {
	Jdet, Jinv        []utils.Matrix    // Sharded mesh Jacobian and inverse transform
	Q_Face            [][4]utils.Matrix // Solution values stored at edge points of RT element
	RHSQ, Residual    [][4]utils.Matrix // Solution Residual storage
	F_RT_DOF          [][4]utils.Matrix // Scalar (projected) flux used for divergence, on all RT element points
	FluxJacInv        [][]utils.Matrix  // Flux Jacobian inverse, one 4x4 matrix (2D) per int point of RT element
	GlobalDT, Time    float64
	Kmax              []int // Local element count (dimension: Kmax[ParallelDegree])
	Np, Nedge, NpFlux int   // Number of points in solution, edge and flux total
	toWorker          []chan struct{}
	fromWorkers       chan int8
}

func (c *Euler) NewElementImplicit() (ei *ElementImplicit) {
	var (
		pm   = c.Partitions
		NPar = pm.ParallelDegree
	)
	ei = &ElementImplicit{
		Jdet:        c.ShardByKTranspose(c.dfr.Jdet),
		Jinv:        c.ShardByKTranspose(c.dfr.Jinv),
		Q_Face:      make([][4]utils.Matrix, NPar),
		RHSQ:        make([][4]utils.Matrix, NPar),
		Residual:    make([][4]utils.Matrix, NPar),
		F_RT_DOF:    make([][4]utils.Matrix, NPar),
		FluxJacInv:  make([][]utils.Matrix, NPar),
		Kmax:        make([]int, NPar),
		Np:          c.dfr.SolutionElement.Np,
		Nedge:       c.dfr.FluxElement.Nedge,
		NpFlux:      c.dfr.FluxElement.Np,
		fromWorkers: make(chan int8, NPar),
		toWorker:    make([]chan struct{}, NPar),
	}
	for np := 0; np < NPar; np++ {
		ei.toWorker[np] = make(chan struct{}, 1)
	}
	for np := 0; np < NPar; np++ {
		// One flux inverse matrix stored per solution point
		ei.FluxJacInv[np] = make([]utils.Matrix, ei.Np*ei.Kmax[np])
	}
	// Initialize memory
	for np := 0; np < NPar; np++ {
		ei.Kmax[np] = pm.GetBucketDimension(np)
		for n := 0; n < 4; n++ {
			ei.Residual[np][n] = utils.NewMatrix(ei.Np, ei.Kmax[np])
			ei.RHSQ[np][n] = utils.NewMatrix(ei.Np, ei.Kmax[np])
			ei.F_RT_DOF[np][n] = utils.NewMatrix(ei.NpFlux, ei.Kmax[np])
			ei.Q_Face[np][n] = utils.NewMatrix(ei.Nedge*3, ei.Kmax[np])
		}
		for i := 0; i < ei.Np*ei.Kmax[np]; i++ { // One 4x4 (2D) flux inverse matrix stored per solution point
			ei.FluxJacInv[np][i] = utils.NewMatrix(4, 4)
		}
	}
	return
}

func (ei *ElementImplicit) StartWorkers(c *Euler) {
	var (
		pm = c.Partitions
		NP = pm.ParallelDegree
	)
	for np := 0; np < NP; np++ {
		go ei.StepWorker(c, np, ei.toWorker[np], ei.fromWorkers)
	}
}

func (ei *ElementImplicit) Step(c *Euler) {
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
			ei.toWorker[np] <- struct{}{}
		}
		currentStep := int8(-1)
		var ts int8
		for np := 0; np < NP; np++ {
			// Each worker sends it's completed step number, we check to make sure they are all the same
			ts = <-ei.fromWorkers
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
		case currentStep < 0:
			return
		}
	}
}

func (ei *ElementImplicit) WorkerDone(subStepP *int8, fromWorker chan int8, isDone bool) {
	if isDone {
		*subStepP = -1
	} else {
		*subStepP++
	}
	fromWorker <- *subStepP // Inform controller what step we've just completed
}

func (ei *ElementImplicit) StepWorker(c *Euler, myThread int, fromController chan struct{}, toController chan int8) {
	// Implements an Element Jacobi time advancement for steady state problems
	var (
		Q0                         = c.Q
		Np                         = ei.Np
		Kmax, Jdet, Jinv, F_RT_DOF = ei.Kmax, ei.Jdet, ei.Jinv, ei.F_RT_DOF
		Q_Face, RHSQ               = ei.Q_Face, ei.RHSQ
		dT                         float64
		subStep                    int8
	)
	for {
		_ = <-fromController // Block until parent sends "go"
		c.PrepareEdgeFlux(Kmax[myThread], Jdet[myThread], Jinv[myThread], F_RT_DOF[myThread], Q0[myThread], Q_Face[myThread])
		ei.WorkerDone(&subStep, toController, false)

		_ = <-fromController                                                                                      // Block until parent sends "go"
		c.SetNormalFluxOnEdges(ei.Time, false, Jdet, nil, F_RT_DOF, Q_Face, c.SortedEdgeKeys[myThread], nil, nil) // Global
		ei.WorkerDone(&subStep, toController, false)

		_ = <-fromController // Block until parent sends "go"
		c.RHS(Kmax[myThread], Jdet[myThread], F_RT_DOF[myThread], RHSQ[myThread])
		dT = ei.GlobalDT
		for n := 0; n < 4; n++ {
			for i := 0; i < Kmax[myThread]*Np; i++ {
				_ = dT
				//Q1[myThread][n].DataP[i] = Q0[myThread][n].DataP[i] + 0.5*RHSQ[myThread][n].DataP[i]*dT
			}
		}
		ei.WorkerDone(&subStep, toController, true)
	}
}

func (ei *ElementImplicit) SetFluxJacobian(c *Euler, Kmax int, Jdet, Jinv utils.Matrix, F_RT_DOF, Q, Q_Face [4]utils.Matrix) {
	// Calculates the flux jacobian for all interior points within the RT element
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
	//dFluxdU - dFdU+dGdU sum of flux jacobians
	// w1 and w2 are the basis function normal vectors: [w1,w2]
	// 0-Nint: [1,0], Nint-2Nint: [0,1], 2Nint-2Nint+Nedge: [w1edge,w2edge]
	//T[0][1] = w1
	//T[0][2] = w2
	//T[1][0] = -w1*(u*u-(1.0/(rho*rho)*((rho*rho)*(u*u)+(rho*rho)*(v*v))*(gamma-1.0))/2.0) - u*v*w2
	//T[1][1] = v*w2 - u*w1*(gamma-3.0)
	//T[1][2] = u*w2 - v*w1*(gamma-1.0)
	//T[1][3] = w1 * (gamma - 1.0)
	//T[2][0] = -w2*(v*v-(1.0/(rho*rho)*((rho*rho)*(u*u)+(rho*rho)*(v*v))*(gamma-1.0))/2.0) - u*v*w1
	//T[2][1] = v*w1 - u*w2*(gamma-1.0)
	//T[2][2] = u*w1 - v*w2*(gamma-3.0)
	//T[2][3] = w2 * (gamma - 1.0)
	//T[3][0] = -((u*w1 + v*w2) * (E*gamma + rho*(u*u) + rho*(v*v) - gamma*rho*(u*u) - gamma*rho*(v*v))) / rho
	//T[3][1] = w1*((E+(E-((rho*rho)*(u*u)+(rho*rho)*(v*v))/(rho*2.0))*(gamma-1.0))/rho-(u*u)*(gamma-1.0)) - u*v*w2*(gamma-1.0)
	//T[3][2] = w2*((E+(E-((rho*rho)*(u*u)+(rho*rho)*(v*v))/(rho*2.0))*(gamma-1.0))/rho-(v*v)*(gamma-1.0)) - u*v*w1*(gamma-1.0)
	//T[3][3] = gamma * (u*w1 + v*w2)
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
		Q0                           = c.Q[myThread]
		Np                           = rk.Np
		Kmax, Jdet, Jinv, F_RT_DOF   = rk.Kmax[myThread], rk.Jdet[myThread], rk.Jinv[myThread], rk.F_RT_DOF[myThread]
		DT, Q_Face, Q1, Q2, Q3, RHSQ = rk.DT[myThread], rk.Q_Face[myThread], rk.Q1[myThread], rk.Q2[myThread], rk.Q3[myThread], rk.RHSQ[myThread]
		Residual                     = rk.Residual[myThread]
		SortedEdgeKeys               = c.SortedEdgeKeys[myThread]
		dT                           float64
		subStep                      int8
		preCon                       = false
		Nedge                        = c.dfr.FluxElement.Nedge
		EdgeQ1, EdgeQ2               = make([][4]float64, Nedge), make([][4]float64, Nedge) // Local working memory
	)
	if c.LocalTimeStepping {
		preCon = false
	}
	for {
		_ = <-fromController // Block until parent sends "go"
		if c.LocalTimeStepping {
			// Setup local time stepping
			for k := 0; k < Kmax; k++ {
				DT.DataP[k] = -100 // Global
			}
		}
		c.PrepareEdgeFlux(Kmax, Jdet, Jinv, F_RT_DOF, Q0, Q_Face)
		rk.WorkerDone(&subStep, toController, false)

		_ = <-fromController // Block until parent sends "go"
		rk.MaxWaveSpeed[myThread] = c.SetNormalFluxOnEdges(rk.Time, true,
			rk.Jdet, rk.DT, rk.F_RT_DOF, rk.Q_Face, SortedEdgeKeys, EdgeQ1, EdgeQ2) // Global
		rk.WorkerDone(&subStep, toController, false)

		_ = <-fromController // Block until parent sends "go"
		if c.LocalTimeStepping {
			// Replicate local time step to the other solution points for each k
			for k := 0; k < Kmax; k++ {
				DT.DataP[k] = c.CFL / DT.DataP[k] // Set each element's DT to CFL/(max_wave_speed)
			}
			// Set the DT of all interior points of each element to the element DT
			for i := 1; i < c.dfr.SolutionElement.Np; i++ {
				for k := 0; k < Kmax; k++ {
					ind := k + Kmax*i
					DT.DataP[ind] = DT.DataP[k]
				}
			}
		}
		c.RHS(Kmax, Jdet, F_RT_DOF, RHSQ)
		if preCon {
			c.PreconditionRHS(Q0, RHSQ, DT, true)
		}
		dT = rk.GlobalDT
		for n := 0; n < 4; n++ {
			for i := 0; i < Kmax*Np; i++ {
				if c.LocalTimeStepping {
					dT = DT.DataP[i]
				}
				Q1[n].DataP[i] = Q0[n].DataP[i] + 0.5*RHSQ[n].DataP[i]*dT
			}
		}
		c.PrepareEdgeFlux(Kmax, Jdet, Jinv, F_RT_DOF, Q1, Q_Face)
		rk.WorkerDone(&subStep, toController, false)

		_ = <-fromController // Block until parent sends "go"
		c.SetNormalFluxOnEdges(rk.Time, false,
			rk.Jdet, rk.DT, rk.F_RT_DOF, rk.Q_Face, SortedEdgeKeys, EdgeQ1, EdgeQ2) // Global
		rk.WorkerDone(&subStep, toController, false)

		_ = <-fromController // Block until parent sends "go"
		c.RHS(Kmax, Jdet, F_RT_DOF, RHSQ)
		if preCon {
			c.PreconditionRHS(Q1, RHSQ, DT, false)
		}
		dT = rk.GlobalDT
		for n := 0; n < 4; n++ {
			for i := 0; i < Kmax*Np; i++ {
				if c.LocalTimeStepping {
					dT = DT.DataP[i]
				}
				Q2[n].DataP[i] = Q1[n].DataP[i] + 0.25*RHSQ[n].DataP[i]*dT
			}
		}
		c.PrepareEdgeFlux(Kmax, Jdet, Jinv, F_RT_DOF, Q2, Q_Face)
		rk.WorkerDone(&subStep, toController, false)

		_ = <-fromController // Block until parent sends "go"
		c.SetNormalFluxOnEdges(rk.Time, false,
			rk.Jdet, rk.DT, rk.F_RT_DOF, rk.Q_Face, SortedEdgeKeys, EdgeQ1, EdgeQ2) // Global
		rk.WorkerDone(&subStep, toController, false)

		_ = <-fromController // Block until parent sends "go"
		c.RHS(Kmax, Jdet, F_RT_DOF, RHSQ)
		if preCon {
			c.PreconditionRHS(Q2, RHSQ, DT, false)
		}
		dT = rk.GlobalDT
		for n := 0; n < 4; n++ {
			for i := 0; i < Kmax*Np; i++ {
				if c.LocalTimeStepping {
					dT = DT.DataP[i]
				}
				Q3[n].DataP[i] = (1. / 3.) * (2*Q0[n].DataP[i] + Q2[n].DataP[i] + RHSQ[n].DataP[i]*dT)
			}
		}
		c.PrepareEdgeFlux(Kmax, Jdet, Jinv, F_RT_DOF, Q3, Q_Face)
		rk.WorkerDone(&subStep, toController, false)

		_ = <-fromController // Block until parent sends "go"
		c.SetNormalFluxOnEdges(rk.Time, false,
			rk.Jdet, rk.DT, rk.F_RT_DOF, rk.Q_Face, SortedEdgeKeys, EdgeQ1, EdgeQ2) // Global
		rk.WorkerDone(&subStep, toController, false)

		_ = <-fromController // Block until parent sends "go"
		c.RHS(Kmax, Jdet, F_RT_DOF, RHSQ)
		if preCon {
			c.PreconditionRHS(Q3, RHSQ, DT, false)
		}
		// Note, we are re-using Q1 as storage for Residual here
		dT = rk.GlobalDT
		for n := 0; n < 4; n++ {
			for i := 0; i < Kmax*Np; i++ {
				if c.LocalTimeStepping {
					dT = DT.DataP[i]
				}
				Residual[n].DataP[i] = Q3[n].DataP[i] + 0.25*RHSQ[n].DataP[i]*dT - Q0[n].DataP[i]
				Q0[n].DataP[i] += Residual[n].DataP[i]
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
		c.dfr.FluxElement.DivInt.Mul(F_RT_DOF[n], RHSQ[n])
		c.DivideByJacobian(Kmax, c.dfr.FluxElement.Nint, Jdet, RHSQ[n].DataP, -1)
	}
	return
}

func (c *Euler) PreconditionRHS(Q, RHSQ [4]utils.Matrix, DT utils.Matrix, calcDT bool) {
	var (
		newRHS       [4]float64
		P0           [4][4]float64
		maxWave, UPC float64
	)
	for i, rho := range Q[0].DataP {
		rhoU, rhoV, E := Q[1].DataP[i], Q[2].DataP[i], Q[3].DataP[i]
		//P0, maxWave, UPC = c.GetPreconditioner(rho, rhoU, rhoV, E)
		P0, maxWave, UPC = c.GetPreconditionerWithError(rho, rhoU, rhoV, E)
		for rNum := 0; rNum < 4; rNum++ {
			newRHS[rNum] =
				P0[rNum][0]*RHSQ[0].DataP[i] +
					P0[rNum][1]*RHSQ[1].DataP[i] +
					P0[rNum][2]*RHSQ[2].DataP[i] +
					P0[rNum][3]*RHSQ[3].DataP[i]
		}
		for rNum := 0; rNum < 4; rNum++ {
			RHSQ[rNum].DataP[i] = newRHS[rNum]
		}
		if calcDT {
			// Unpack base local time step to obtain jacobian and order metric
			met := c.CFL / DT.DataP[i] / UPC
			DT.DataP[i] = c.CFL / (met * maxWave)
		}
	}
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
func (c *Euler) PrintUpdate(Time, dt float64, steps int, Q, Residual [][4]utils.Matrix, plotQ bool, pm *PlotMeta,
	printMem bool) {
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
	if printMem {
		fmt.Printf(" :: %s", utils.GetMemUsage())
	}
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

func (c *Euler) GetPreconditioner(rho, rhoU, rhoV, E float64) (P0 [4][4]float64, maxWave, UPC float64) {
	/*
		[nx,ny] is the direction vector of the flux at the evaluated point:
			Flux = [F,G] = |[F,G]| * [nx,ny]
	*/
	var (
		gamma  = c.FS.Gamma
		sqrt   = math.Sqrt
		GM1    = gamma - 1.0
		u, v   = rhoU / rho, rhoV / rho
		u2, v2 = u * u, v * v
		qq     = u2 + v2
		uave   = sqrt(qq)
		qqq    = GM1 * qq / 2
		P      = GM1 * (E - 0.5*qq*rho)
		P2     = P * P
		c2     = gamma * P / rho
		C      = sqrt(c2)
		mach   = uave / C
		m2     = mach * mach
		h      = c2/GM1 + qq/2
		K1     = 1.0 // should be between 1 and 1.1
		Minf   = c.FS.Minf
		Minf2  = Minf * Minf
		Minf4  = Minf2 * Minf2
		B2     = K1 * m2 * (1 + (1-K1*Minf2)*m2/(K1*Minf4))
	)
	B2 = math.Max(1, B2)
	maxWave = 0.5 * ((B2+1)*uave + sqrt((B2-1)*(B2-1)*qq+4*B2*c2))
	UPC = uave + C
	P0[0][0] = 1.0/P2 + (1.0/(P2)*qqq*(B2*P2-1.0))/c2
	P0[0][1] = -(1.0 / P2 * u * (B2*P2 - 1.0) * GM1) / c2
	P0[0][2] = -(1.0 / P2 * v * (B2*P2 - 1.0) * GM1) / c2
	P0[0][3] = (1.0 / P2 * (B2*P2 - 1.0) * GM1) / c2
	P0[1][0] = -1.0/P2*u*(P2-1.0) + (1.0/P2*qqq*u*(B2*P2-1.0))/c2
	P0[1][1] = -(1.0/P2*(u*u)*(B2*P2-1.0)*GM1)/c2 + 1.0
	P0[1][2] = -(1.0 / P2 * u * v * (B2*P2 - 1.0) * GM1) / c2
	P0[1][3] = (1.0 / P2 * u * (B2*P2 - 1.0) * GM1) / c2
	P0[2][0] = -1.0/P2*v*(P2-1.0) + (1.0/P2*qqq*v*(B2*P2-1.0))/c2
	P0[2][1] = -(1.0 / P2 * u * v * (B2*P2 - 1.0) * GM1) / c2
	P0[2][2] = -(1.0/P2*(v*v)*(B2*P2-1.0)*GM1)/c2 + 1.0
	P0[2][3] = (1.0 / P2 * v * (B2*P2 - 1.0) * GM1) / c2
	P0[3][0] = -(c2*(u*u)+c2*(v*v)-B2*h*qqq)/c2 + (1.0/P2*m2*(c2-qqq))/2.0
	P0[3][1] = (1.0/P2*u*(-m2+gamma*m2+P2*2.0))/2.0 - (B2*h*u*GM1)/c2
	P0[3][2] = (1.0/P2*v*(-m2+gamma*m2+P2*2.0))/2.0 - (B2*h*v*GM1)/c2
	P0[3][3] = 1.0/P2*m2*GM1*(-1.0/2.0) + (B2*h*GM1)/c2
	return
}

func (c *Euler) GetPreconditionerWithError(rho, rhoU, rhoV, E float64) (P0 [4][4]float64, maxWave, UPC float64) {
	/*
			[nx,ny] is the direction vector of the flux at the evaluated point:
				Flux = [F,G] = |[F,G]| * [nx,ny]

		This version includes an error in the Turkel paper where the last row in the variable transform Jacobian
		from Entropy variables to Conservative variables is divided by P
	*/
	var (
		gamma  = c.FS.Gamma
		sqrt   = math.Sqrt
		GM1    = gamma - 1.0
		u, v   = rhoU / rho, rhoV / rho
		u2, v2 = u * u, v * v
		qq     = u2 + v2
		uave   = sqrt(qq)
		qqq    = GM1 * qq / 2
		P      = GM1 * (E - 0.5*qq*rho)
		c2     = gamma * P / rho
		C      = sqrt(c2)
		mach   = uave / C
		m2     = mach * mach
		h      = c2/GM1 + qq/2
		K1     = 1.0 // should be between 1 and 1.1
		Minf   = c.FS.Minf
		Minf2  = Minf * Minf
		Minf4  = Minf2 * Minf2
		B2     = K1 * m2 * (1 + (1-K1*Minf2)*m2/(K1*Minf4))
	)
	B2 = math.Max(1, B2)
	maxWave = 0.5 * ((B2+1)*uave + sqrt((B2-1)*(B2-1)*qq+4*B2*c2))
	UPC = uave + C
	P0[0][0] = 1.0/P + (qqq*(B2*P-1.0))/(P*c2)
	P0[0][1] = -(u * (B2*P - 1.0) * (gamma - 1.0)) / (P * c2)
	P0[0][2] = -(v * (B2*P - 1.0) * (gamma - 1.0)) / (P * c2)
	P0[0][3] = ((B2*P - 1.0) * (gamma - 1.0)) / (P * c2)
	P0[1][0] = -(u*(P-1.0))/P + (qqq*u*(B2*P-1.0))/(P*c2)
	P0[1][1] = -((u*u)*(B2*P-1.0)*(gamma-1.0))/(P*c2) + 1.0
	P0[1][2] = -(u * v * (B2*P - 1.0) * (gamma - 1.0)) / (P * c2)
	P0[1][3] = (u * (B2*P - 1.0) * (gamma - 1.0)) / (P * c2)
	P0[2][0] = -(v*(P-1.0))/P + (qqq*v*(B2*P-1.0))/(P*c2)
	P0[2][1] = -(u * v * (B2*P - 1.0) * (gamma - 1.0)) / (P * c2)
	P0[2][2] = -((v*v)*(B2*P-1.0)*(gamma-1.0))/(P*c2) + 1.0
	P0[2][3] = (v * (B2*P - 1.0) * (gamma - 1.0)) / (P * c2)
	P0[3][0] = -u*u - v*v + (m2*(c2-qqq))/(P*2.0) + (B2*h*qqq)/c2
	P0[3][1] = (u*(P*2.0-m2+gamma*m2))/(P*2.0) - (B2*h*u*(gamma-1.0))/c2
	P0[3][2] = (v*(P*2.0-m2+gamma*m2))/(P*2.0) - (B2*h*v*(gamma-1.0))/c2
	P0[3][3] = (m2*(gamma-1.0)*(-1.0/2.0))/P + (B2*h*(gamma-1.0))/c2
	return
}

package Euler2D

import (
	"fmt"
	"math"
	"sync"
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

	c.PartitionEdges() // Setup the key for edge calculations, useful for parallelizing the process

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
		Time, dt  float64
		steps     int
		finished  bool
		plotQ     = pm.Plot
		Residual  [][4]utils.Matrix
	)
	if c.profile {
		//defer profile.Start(profile.CPUProfile).Stop()
		defer profile.Start().Stop()
	}

	c.PrintInitialization(FinalTime)

	rk := c.NewRungeKuttaSSP()

	elapsed := time.Duration(0)
	var start time.Time
	for !finished {
		start = time.Now()
		Residual, dt = rk.Step(c, Time, c.Q)
		elapsed += time.Now().Sub(start)
		steps++
		Time += dt
		finished = c.CheckIfFinished(Time, FinalTime, steps)
		if finished || steps%pm.StepsBeforePlot == 0 || steps == 1 {
			c.PrintUpdate(Time, dt, steps, c.Q, Residual, plotQ, pm)
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
	Kmax              []int             // Local element count (dimension: Kmax[ParallelDegree])
	Np, Nedge, NpFlux int               // Number of points in solution, edge and flux total
}

func (c *Euler) NewRungeKuttaSSP() (rk *RungeKutta4SSP) {
	var (
		pm   = c.Partitions
		NPar = pm.ParallelDegree
	)
	rk = &RungeKutta4SSP{
		Jdet:     c.ShardByKTranspose(c.dfr.Jdet),
		Jinv:     c.ShardByKTranspose(c.dfr.Jinv),
		RHSQ:     make([][4]utils.Matrix, NPar),
		Q_Face:   make([][4]utils.Matrix, NPar),
		Q1:       make([][4]utils.Matrix, NPar),
		Q2:       make([][4]utils.Matrix, NPar),
		Q3:       make([][4]utils.Matrix, NPar),
		Residual: make([][4]utils.Matrix, NPar),
		F_RT_DOF: make([][4]utils.Matrix, NPar),
		DT:       make([]utils.Matrix, NPar),
		Kmax:     make([]int, NPar),
		Np:       c.dfr.SolutionElement.Np,
		Nedge:    c.dfr.FluxElement.Nedge,
		NpFlux:   c.dfr.FluxElement.Np,
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
	return
}

func (rk *RungeKutta4SSP) Step(c *Euler, Time float64, Q0 [][4]utils.Matrix) (Residual [][4]utils.Matrix, dt float64) {
	var (
		Np                           = rk.Np
		Kmax, Jdet, Jinv, F_RT_DOF   = rk.Kmax, rk.Jdet, rk.Jinv, rk.F_RT_DOF
		DT, Q_Face, Q1, Q2, Q3, RHSQ = rk.DT, rk.Q_Face, rk.Q1, rk.Q2, rk.Q3, rk.RHSQ
		pm                           = c.Partitions
		NP                           = pm.ParallelDegree
		wg                           = sync.WaitGroup{}
	)
	Residual = Q1
	for np := 0; np < NP; np++ {
		wg.Add(1)
		go func(np int) {
			c.PrepareEdgeFlux(Kmax[np], Jdet[np], Jinv[np], F_RT_DOF[np], Q0[np], Q_Face[np])
			wg.Done()
		}(np)
	}
	wg.Wait()
	dt = c.ParallelEdgeUpdate(Time, true, Jdet, DT, F_RT_DOF, Q_Face) // Must sync parallel before calling
	if Time+dt > c.FinalTime {
		dt = c.FinalTime - Time
	}
	for np := 0; np < NP; np++ {
		wg.Add(1)
		go func(np int) {
			c.RHS(Kmax[np], Jdet[np], F_RT_DOF[np], RHSQ[np])
			dT := dt
			for n := 0; n < 4; n++ {
				for i := 0; i < Kmax[np]*Np; i++ {
					if c.LocalTimeStepping {
						dT = DT[np].DataP[i]
					}
					Q1[np][n].DataP[i] = Q0[np][n].DataP[i] + 0.5*RHSQ[np][n].DataP[i]*dT
				}
			}
			c.PrepareEdgeFlux(Kmax[np], Jdet[np], Jinv[np], F_RT_DOF[np], Q1[np], Q_Face[np])
			wg.Done()
		}(np)
	}
	wg.Wait()
	_ = c.ParallelEdgeUpdate(Time, false, Jdet, DT, F_RT_DOF, Q_Face) // Must sync parallel before calling
	for np := 0; np < NP; np++ {
		wg.Add(1)
		go func(np int) {
			c.RHS(Kmax[np], Jdet[np], F_RT_DOF[np], RHSQ[np])
			dT := dt
			for n := 0; n < 4; n++ {
				for i := 0; i < Kmax[np]*Np; i++ {
					if c.LocalTimeStepping {
						dT = DT[np].DataP[i]
					}
					Q2[np][n].DataP[i] = Q1[np][n].DataP[i] + 0.25*RHSQ[np][n].DataP[i]*dT
				}
			}
			c.PrepareEdgeFlux(Kmax[np], Jdet[np], Jinv[np], F_RT_DOF[np], Q2[np], Q_Face[np])
			wg.Done()
		}(np)
	}
	wg.Wait()
	_ = c.ParallelEdgeUpdate(Time, false, Jdet, DT, F_RT_DOF, Q_Face) // Must sync parallel before calling
	for np := 0; np < NP; np++ {
		wg.Add(1)
		go func(np int) {
			c.RHS(Kmax[np], Jdet[np], F_RT_DOF[np], RHSQ[np])
			dT := dt
			for n := 0; n < 4; n++ {
				for i := 0; i < Kmax[np]*Np; i++ {
					if c.LocalTimeStepping {
						dT = DT[np].DataP[i]
					}
					Q3[np][n].DataP[i] = (1. / 3.) * (2*Q0[np][n].DataP[i] + Q2[np][n].DataP[i] + RHSQ[np][n].DataP[i]*dT)
				}
			}
			c.PrepareEdgeFlux(Kmax[np], Jdet[np], Jinv[np], F_RT_DOF[np], Q3[np], Q_Face[np])
			wg.Done()
		}(np)
	}
	wg.Wait()
	_ = c.ParallelEdgeUpdate(Time, false, Jdet, DT, F_RT_DOF, Q_Face) // Must sync parallel before calling
	for np := 0; np < NP; np++ {
		wg.Add(1)
		go func(np int) {
			c.RHS(Kmax[np], Jdet[np], F_RT_DOF[np], RHSQ[np])
			// Note, we are re-using Q1 as storage for Residual here
			dT := dt
			for n := 0; n < 4; n++ {
				for i := 0; i < Kmax[np]*Np; i++ {
					if c.LocalTimeStepping {
						dT = DT[np].DataP[i]
					}
					Residual[np][n].DataP[i] = Q3[np][n].DataP[i] + 0.25*RHSQ[np][n].DataP[i]*dT - Q0[np][n].DataP[i]
					Q0[np][n].DataP[i] += Residual[np][n].DataP[i]
				}
			}
			wg.Done()
		}(np)
	}
	wg.Wait()
	return
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
}

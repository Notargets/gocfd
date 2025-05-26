package Euler2D

import (
	"fmt"
	"math"
	"os"
	"runtime"
	"runtime/pprof"
	"sync"
	"syscall"
	"time"

	"github.com/notargets/gocfd/InputParameters"

	"github.com/notargets/gocfd/model_problems/Euler2D/sod_shock_tube"

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
	MeshFile           string
	CFL, FinalTime     float64
	FSFar, FSIn, FSOut *FreeStream
	DFR                *DG2D.DFR2D
	profile            bool // Generate a CPU profile of the solver
	FluxCalcAlgo       FluxType
	Case               InitType
	AnalyticSolution   ExactState
	FluxCalcMock       func(rho, rhoU, rhoV, E float64) (Fx, Fy [4]float64) // For testing
	SortedEdgeKeys     []EdgeKeySlice                                       // Buckets, one for each parallel partition
	Partitions         *utils.PartitionMap                                  // mapping of elements into bins for parallelism
	LocalTimeStepping  bool
	MaxIterations      int
	// Below are partitioned by K (elements) in the first slice
	Q                    [][4]utils.Matrix // Sharded solution variables, stored at solution point locations, Np_solution x K
	Q4                   [4]utils.Matrix   // Non-Sharded solution for plotting
	SolutionX, SolutionY []utils.Matrix
	ShockFinder          *DG2D.ModeAliasShockFinder
	Dissipation          *ScalarDissipation
	Kappa                float64 // Artificial Dissipation Strength constant
	// Edge number mapped quantities, i.e. Face Normal Flux
	EdgeStore           *EdgeValueStorage
	ShockTube           *sod_shock_tube.SODShockTube
	SolutionFieldWriter *DG2D.AVSFieldWriter
	RK                  *RungeKutta5SSP
	NeighborNotifier    *utils.NeighborNotifier
}

func NewEuler(ip *InputParameters.InputParameters2D, meshFile string, ProcLimit int, verbose, profile bool) (c *Euler) {
	c = &Euler{
		MeshFile:          meshFile,
		CFL:               ip.CFL,
		FinalTime:         ip.FinalTime,
		FluxCalcAlgo:      NewFluxType(ip.FluxType),
		Case:              NewInitType(ip.InitType),
		LocalTimeStepping: ip.LocalTimeStepping,
		MaxIterations:     ip.MaxIterations,
		FSFar:             NewFreeStream(ip.Minf, ip.Gamma, ip.Alpha),
		profile:           profile,
		Kappa:             ip.Kappa,
	}
	c.FluxCalcMock = c.FluxCalcBase

	if len(meshFile) == 0 {
		return
	}

	// Read mesh file, initialize geometry and finite elements
	c.DFR = DG2D.NewDFR2D(ip.PolynomialOrder, verbose, meshFile)
	c.ShockFinder = c.DFR.NewAliasShockFinder(c.Kappa)
	fmd := &DG2D.FieldMetadata{
		NumFields:  1,
		FieldNames: ip.PlotFields,
		SolutionMetadata: map[string]interface{}{
			"CFL":     ip.CFL,
			"Mesh":    meshFile,
			"Case":    c.Case.String(),
			"Solver":  "2D Euler",
			"FS Mach": ip.Minf,
		},
		GitVersion: "",
	}
	c.SolutionFieldWriter = c.DFR.NewAVSFieldWriter(
		fmd, "solutionfile.gobcfd", c.DFR.GraphMesh)

	c.SetParallelDegree(ProcLimit, c.DFR.K) // Must occur after determining the number of elements
	c.PartitionEdgesByK()                   // Setup the key for edge calculations, useful for parallelizing the process
	c.NeighborNotifier = utils.NewNeighborNotifier(c.Partitions, c.DFR.Tris.EtoE)

	// Allocate Normal flux storage and indices
	c.EdgeStore = c.NewEdgeStorage()

	c.InitializeSolution(verbose)

	// Allocate a solution limiter
	lt := NewLimiterType(ip.Limiter)
	// Initiate Artificial Dissipation
	if lt == PerssonC0T && c.DFR.N != 0 {
		c.Dissipation = NewScalarDissipation(c.Kappa, C0, c.DFR, c.Partitions)
	}

	// Save graph mesh
	c.DFR.OutputMesh("meshfile.gobcfd", c.SerializeBCs())

	c.RK = c.NewRungeKuttaSSP()

	if verbose {
		fmt.Printf("Euler Equations in 2 Dimensions\n")
		fmt.Printf("Using %d go routines in parallel\n", c.Partitions.ParallelDegree)
		fmt.Printf("Solving %s\n", c.Case.String())
		switch c.Case {
		case FREESTREAM:
			fmt.Printf("Mach Infinity = %8.5f, Angle of Attack = %8.5f\n", ip.Minf, ip.Alpha)
		}
		fmt.Printf("Flux Algorithm: [%s]\n", c.FluxCalcAlgo.Print())
		if c.Dissipation != nil {
			fmt.Printf("Artificial Dissipation Coefficient: Kappa = [%5.3f]\n", c.Dissipation.Kappa)
		}
		fmt.Printf("Shock Finder Coefficient: Kappa = [%5.3f]\n", c.Kappa)
		fmt.Printf("CFL = %8.4f, Polynomial Degree N = %d (1 is linear), Num Elements K = %d\n\n\n",
			ip.CFL, ip.PolynomialOrder, c.DFR.K)
	}
	return
}

func (c *Euler) Solve() {
	var (
		FinalTime = c.FinalTime
		steps     int
		finished  bool
	)
	if c.profile {
		// defer profile.Start().Stop()
		// Start CPU profile
		fCPU, err := os.Create("cpu.pprof")
		if err != nil {
			panic(fmt.Errorf("could not create CPU profile: %s", err))
		}
		pprof.StartCPUProfile(fCPU)
		defer func() {
			pprof.StopCPUProfile()
			fCPU.Close()
		}()

		// Run your workload here

		// After workload, write heap profile
		fMem, err := os.Create("mem.pprof")
		if err != nil {
			panic(fmt.Errorf("could not create memory profile: %s", err))
		}
		runtime.GC() // run GC to get up-to-date mem profile
		if err := pprof.WriteHeapProfile(fMem); err != nil {
			panic(fmt.Errorf("could not write memory profile: %s", err))
		}
		fMem.Close()
	}

	c.PrintInitialization(FinalTime)

	elapsed := time.Duration(0)
	var start time.Time
	for !finished {
		start = time.Now()
		c.RK.Step(c)
		elapsed += time.Now().Sub(start)
		steps++
		c.RK.Time += c.RK.GlobalDT
		c.RK.StepCount++
		finished = c.CheckIfFinished(c.RK.Time, steps)
		// if finished || steps == 1 || steps%1 == 0 {
		if finished || steps == 1 || steps%100 == 0 {
			var printMem bool
			if steps%100 == 0 {
				printMem = true
			}
			c.PrintUpdate(c.RK.Time, c.RK.GlobalDT, steps, c.Q, c.RK.Residual, printMem, c.RK.LimitedPoints)
			if len(c.SolutionFieldWriter.FieldMeta.FieldNames) != 0 {
				c.RecombineShardsKBy4(c.Q, &c.Q4)
				fieldMap := make(map[string][]float64)
				var nFields, length int
				for _, name := range c.SolutionFieldWriter.FieldMeta.FieldNames {
					ff, match := BestMatchFlowFunction(name)
					if !match {
						panic("Unable to find matching flow function named: " + name)
						// } else {
						// 	fmt.Printf("Field function: %s\n", ff.String())
					}
					field := c.GetPlotField(c.Q4, ff)
					fieldMap[name] = field.DataP
					length = len(fieldMap[name])
					nFields++
					// fmt.Printf("Writing field:%s, Min/Max=%.2f/%.2f\n",
					// 	name, FMat.Min(), FMat.Max())
				}
				c.SolutionFieldWriter.Save(&DG2D.SingleFieldMetadata{
					Iterations: steps,
					Time:       float32(c.RK.Time),
					Count:      nFields,
					Length:     length,
				},
					fieldMap)
			}
		}
	}
	c.PrintFinal(elapsed, steps)
	c.OutputFinal()
}

func (c *Euler) OutputFinal() {
	switch c.Case {
	case SHOCKTUBE:
		var (
			st = c.ShockTube
		)
		c.RecombineShardsKBy4(c.Q, &c.Q4)
		c.ShockTube.InterpolateFields(c.Q4)
		file, err := os.Create("shocktube.dat")
		if err != nil {
			panic(fmt.Errorf("could not create shocktube.dat: %s", err))
		}
		file.WriteString(fmt.Sprintf("Meshfile: %s\n", c.MeshFile))
		file.WriteString(fmt.Sprintf("X\tRho\tRhoU\tE\n"))
		for i, x := range st.XLocations {
			file.WriteString(fmt.Sprintf("%.8f\t%.8f\t%.8f\t%.8f\n",
				x, st.Rho[i], st.RhoU[i], st.E[i]))
		}
		file, err = os.Create("shocktube_analytic.dat")
		if err != nil {
			panic(fmt.Errorf("could not create shocktube_analytic.dat: %s", err))
		}
		file.WriteString(fmt.Sprintf("X\tRho\tRhoU\tE\n"))
		XA, rhoA, _, rhoUA, rhoEA :=
			c.ShockTube.GetAnalyticSolution(c.FinalTime)
		for i, x := range XA {
			file.WriteString(fmt.Sprintf("%.8f\t%.8f\t%.8f\t%.8f\n",
				x, rhoA[i], rhoUA[i], rhoEA[i]))
		}
	}
	return
}

type RungeKutta5SSP struct {
	Jdet, Jinv                                 []utils.Matrix       // Sharded mesh Jacobian and inverse transform
	RHSQ, Q_Face                               [][4]utils.Matrix    // State used for matrix multiplies within the time step algorithm
	Flux_Face                                  [][2][4]utils.Matrix // Flux interpolated to edges from interior
	Q1, Q2, Q3, Q4                             [][4]utils.Matrix    // Intermediate solution state
	Residual                                   [][4]utils.Matrix    // Used for reporting, aliased to Q1
	Se                                         []utils.Vector       // Element moment - troubled element indicator
	LScratch                                   [][3]utils.Matrix    // Scratch space for limiter
	F_RT_DOF                                   [][4]utils.Matrix    // Normal flux used for divergence
	DT, DTVisc                                 []utils.Vector       // Local time step storage
	Epsilon                                    []utils.Matrix
	GlobalMaxWaveSpeed, GlobalMaxViscWaveSpeed []float64 // Shard max wavespeed
	GlobalDT, Time                             float64
	StepCount                                  int
	Kmax                                       []int          // Local element count (dimension: Kmax[ParallelDegree])
	NpInt, NpEdge, NpFlux                      int            // Number of points in solution, edge and flux total
	EdgeQ1, EdgeQ2                             [][][4]float64 // Sharded local working memory, dimensions NpEdge
	LimitedPoints                              []int          // Sharded number of limited points
	ShockSensor                                []*DG2D.ModeAliasShockFinder
	VtoE                                       []VertexToElement // Sharded vertex to element map, [2] is [vertID, ElementID_Sharded]
	EtoV                                       []utils.Matrix    // Sharded Element to Vertex map, Kx3
}

func (c *Euler) NewRungeKuttaSSP() (rk *RungeKutta5SSP) {
	var (
		pm   = c.Partitions
		NPar = pm.ParallelDegree
	)
	rk = &RungeKutta5SSP{
		Jdet:                   c.ShardByKTranspose(c.DFR.Jdet),
		Jinv:                   c.ShardByKTranspose(c.DFR.Jinv),
		RHSQ:                   make([][4]utils.Matrix, NPar),
		Q_Face:                 make([][4]utils.Matrix, NPar),
		Flux_Face:              make([][2][4]utils.Matrix, NPar),
		Q1:                     make([][4]utils.Matrix, NPar),
		Q2:                     make([][4]utils.Matrix, NPar),
		Q3:                     make([][4]utils.Matrix, NPar),
		Q4:                     make([][4]utils.Matrix, NPar),
		Residual:               make([][4]utils.Matrix, NPar),
		Se:                     make([]utils.Vector, NPar),
		LScratch:               make([][3]utils.Matrix, NPar),
		F_RT_DOF:               make([][4]utils.Matrix, NPar),
		DT:                     make([]utils.Vector, NPar),
		DTVisc:                 make([]utils.Vector, NPar),
		GlobalMaxWaveSpeed:     make([]float64, NPar),
		GlobalMaxViscWaveSpeed: make([]float64, NPar),
		Kmax:                   make([]int, NPar),
		NpInt:                  c.DFR.SolutionElement.Np,
		NpEdge:                 c.DFR.FluxElement.NpEdge,
		NpFlux:                 c.DFR.FluxElement.Np,
		EdgeQ1:                 make([][][4]float64, NPar),
		EdgeQ2:                 make([][][4]float64, NPar),
		LimitedPoints:          make([]int, NPar),
		ShockSensor:            make([]*DG2D.ModeAliasShockFinder, NPar),
		VtoE:                   NewVertexToElement(c.DFR.Tris.EToV).Shard(pm),
		EtoV:                   c.shardEtoV(c.DFR.Tris.EToV),
	}
	// Initialize memory for RHS
	for np := 0; np < NPar; np++ {
		rk.Kmax[np] = pm.GetBucketDimension(np)
		rk.LScratch[np][0] = utils.NewMatrix(rk.NpInt, rk.Kmax[np])
		rk.LScratch[np][1] = utils.NewMatrix(rk.NpInt, rk.Kmax[np])
		rk.LScratch[np][2] = utils.NewMatrix(rk.NpInt, rk.Kmax[np])
		for n := 0; n < 4; n++ {
			rk.Q1[np][n] = utils.NewMatrix(rk.NpInt, rk.Kmax[np])
			rk.Q2[np][n] = utils.NewMatrix(rk.NpInt, rk.Kmax[np])
			rk.Q3[np][n] = utils.NewMatrix(rk.NpInt, rk.Kmax[np])
			rk.Q4[np][n] = utils.NewMatrix(rk.NpInt, rk.Kmax[np])
			rk.RHSQ[np][n] = utils.NewMatrix(rk.NpInt, rk.Kmax[np])
			rk.Residual[np][n] = utils.NewMatrix(rk.NpInt, rk.Kmax[np])
			rk.F_RT_DOF[np][n] = utils.NewMatrix(rk.NpFlux, rk.Kmax[np])
			rk.Q_Face[np][n] = utils.NewMatrix(rk.NpEdge*3, rk.Kmax[np])
			rk.Flux_Face[np][0][n] = utils.NewMatrix(rk.NpEdge*3, rk.Kmax[np])
			rk.Flux_Face[np][1][n] = utils.NewMatrix(rk.NpEdge*3, rk.Kmax[np])
		}
		rk.Se[np] = utils.NewVector(rk.Kmax[np])
		rk.DT[np] = utils.NewVector(rk.Kmax[np])
		rk.DTVisc[np] = utils.NewVector(rk.Kmax[np])
		rk.EdgeQ1[np] = make([][4]float64, rk.NpEdge)
		rk.EdgeQ2[np] = make([][4]float64, rk.NpEdge)
		rk.ShockSensor[np] = c.DFR.NewAliasShockFinder(c.Kappa)
	}
	if c.Dissipation != nil {
		rk.Epsilon = c.Dissipation.Epsilon
	}
	return
}

func (rk *RungeKutta5SSP) Step(c *Euler) {
	/*
		This is the controller thread - it manages and synchronizes the workers
		NOTE: Any CPU work done in this thread makes the workers wait - be careful!
	*/
	for rkStep := 0; rkStep < 5; rkStep++ {
		// Workers are blocked below here until the StepWorker section - make sure significant work done here is abs necessary!
		// Advance the workers one step
		rk.StepWorker(c, rkStep)
	}
}

func (rk *RungeKutta5SSP) StepWorker(c *Euler, rkStep int) {
	var (
		QQQAll = [][][4]utils.Matrix{c.Q, rk.Q1, rk.Q2, rk.Q3, rk.Q4}[rkStep]
		sum    = func(a, b float64) float64 {
			return a + b
		}
		div3 = func(a float64) float64 {
			return a / 3.
		}
		Max = func(a, b float64) float64 {
			return max(a, b)
		}
	)
	/*
		Inline functions
	*/
	doParallel := func(myFunc func(np int)) {
		wg := sync.WaitGroup{}
		for np := 0; np < c.Partitions.ParallelDegree; np++ {
			wg.Add(1)
			go func(np int, myFunc func(np int)) {
				defer wg.Done()
				myFunc(np)
			}(np, myFunc)
		}
		wg.Wait()
	}
	// For debugging
	doSerial := func(myFunc func(np int)) {
		wg := sync.WaitGroup{}
		for np := 0; np < c.Partitions.ParallelDegree; np++ {
			wg.Add(1)
			func(np int, myFunc func(np int)) {
				defer wg.Done()
				myFunc(np)
			}(np, myFunc)
		}
		wg.Wait()
	}
	_, _ = doSerial, doParallel
	rkAdvance := func(np int) {
		var (
			Kmax, Jdet, Jinv = rk.Kmax[np], rk.Jdet[np], rk.Jinv[np]
			DT               = rk.DT[np]
			F_RT_DOF         = rk.F_RT_DOF[np]
			RHSQ, Residual   = rk.RHSQ[np], rk.Residual[np]
			QQQ              = QQQAll[np]
			// Debug              = false
			Debug = true
			Np, _ = c.Q[np][0].Dims()
		)
		checkFunc := func() {
			if !Debug {
				return
			}
			utils.IsNan(QQQ)
			utils.IsNan(F_RT_DOF)
			utils.IsNan(RHSQ)
			utils.IsNan(Jinv)
			utils.IsNan(Jdet)
			utils.IsNan(c.Q[np])
			utils.IsNan(rk.Q1[np])
			utils.IsNan(rk.Q2[np])
			utils.IsNan(rk.Q3[np])
			utils.IsNan(rk.Q4[np])
			utils.IsNan(DT)
		}
		/*
			SSP54 RK Coefficients from:
			"A numerical study of diagonally split Runge-Kutta methods for PDEs with discontinuities"
			Colin B. Macdonald, Sigal Gottlieb and Steven J. Ruuth, 2007
		*/
		checkFunc()
		c.SetRTFluxInternal(Kmax, Jdet, Jinv, F_RT_DOF, QQQ) // Updates F_RT_DOF with values from Q
		c.SetRTFluxOnEdges(np, Kmax, F_RT_DOF)
		c.RHSInternalPoints(Kmax, Jdet, F_RT_DOF, RHSQ)
		sd := c.Dissipation
		if sd != nil {
			sd.AddDissipation(c, np, Jinv, Jdet, RHSQ)
			sd.LimitFilterSolution(np, RHSQ, rk.LScratch[np][2], rk.ShockSensor[np])
			// sd.FilterSolution(np, RHSQ, rk.LScratch[np][2], rk.ShockSensor[np])
		}
		for n := 0; n < 4; n++ {
			var (
				q0     = c.Q[np][n].DataP
				q1     = rk.Q1[np][n].DataP
				q2     = rk.Q2[np][n].DataP
				q3     = rk.Q3[np][n].DataP
				q4     = rk.Q4[np][n].DataP
				R, RHS = Residual[n].DataP, RHSQ[n].DataP
			)
			switch rkStep {
			case 0:
				for i := 0; i < Np; i++ {
					for k := 0; k < Kmax; k++ {
						ind := k + i*Kmax
						dtRHS := DT.DataP[k] * RHSQ[n].DataP[ind]
						q1[ind] = q0[ind] + 0.391752226571890*dtRHS
					}
				}
			case 1:
				for i := 0; i < Np; i++ {
					for k := 0; k < Kmax; k++ {
						ind := k + i*Kmax
						dtRHS := DT.DataP[k] * RHSQ[n].DataP[ind]
						q2[ind] = 0.444370493651235*q0[ind] +
							0.555629506348765*q1[ind] +
							0.368410593050371*dtRHS
					}
				}
			case 2:
				for i := 0; i < Np; i++ {
					for k := 0; k < Kmax; k++ {
						ind := k + i*Kmax
						dtRHS := DT.DataP[k] * RHSQ[n].DataP[ind]
						q3[ind] = 0.620101851488403*q0[ind] +
							0.379898148511597*q2[ind] +
							0.251891774271694*dtRHS
					}
				}
			case 3:
				for i := 0; i < Np; i++ {
					for k := 0; k < Kmax; k++ {
						ind := k + i*Kmax
						dtRHS := DT.DataP[k] * RHSQ[n].DataP[ind]
						R[ind] = RHS[ind] // Store the current RHS for use in the last RK step
						q4[ind] = 0.178079954393132*q0[ind] +
							0.821920045606868*q3[ind] +
							0.544974750228521*dtRHS
					}
				}
			case 4:
				for i := 0; i < Np; i++ {
					for k := 0; k < Kmax; k++ {
						ind := k + i*Kmax
						dtRHS := DT.DataP[k] * RHSQ[n].DataP[ind]
						dtR3 := DT.DataP[k] * R[ind]
						R[ind] = -q0[ind] + 0.517231671970585*q2[ind] +
							0.096059710526146*q3[ind] +
							0.386708617503269*q4[ind] +
							0.063692468666290*dtR3 + 0.226007483236906*dtRHS
						q0[ind] += R[ind]
					}
				}
			}
		}
	}
	/*
		Time Step Execution
	*/
	var (
		sd = c.Dissipation
	)
	// doSerial(func(np int) {
	doParallel(func(np int) {
		var (
			QQQ = QQQAll[np]
		)
		if sd != nil {
			rk.ShockSensor[np].UpdateSeMoment(QQQ[0], rk.LScratch[np], rk.Se[np])
			sd.UpdateShockFinderSigma(np, rk.Se[np])
			sd.CalculateElementViscosity(np)
		}
	})
	doParallel(func(np int) {
		if sd != nil {
			c.MergeElementScalarToVertices(rk.VtoE[np], sd.SigmaScalar,
				sd.SigmaVertex, Max)
			c.MergeElementScalarToVertices(rk.VtoE[np], sd.EpsilonScalar,
				sd.EpsVertex, Max)
		}
	})
	// doSerial(func(np int) {
	doParallel(func(np int) {
		var (
			QQQ = QQQAll[np]
		)
		if sd != nil {
			c.MergeVertexScalarToElement(sd.SigmaVertex, sd.SigmaScalar[np],
				rk.EtoV[np], sum, div3)

			sd.InterpolateEpsilonSigma(np, rk.EtoV[np])
		}
		// if c.Dissipation != nil {
		// if (rkStep == 1) && c.Dissipation != nil {
		if (rkStep == 2) && c.Dissipation != nil {
			// if (rkStep == 4) && c.Dissipation != nil {
			// if (rkStep > 1) && c.Dissipation != nil {
			sd.LimitFilterSolution(np, QQQ, rk.LScratch[np][2],
				rk.ShockSensor[np], false)
			// sd.LimitSolution(np, QQQ, rk.LScratch[np][2])
		}
		c.InterpolateSolutionToEdges(QQQ, rk.Q_Face[np])
		// if c.Dissipation != nil {
		// 	sd.FilterSolution(np, QQQ, rk.LScratch[np][2], rk.ShockSensor[np])
		// }
	})
	// doSerial(func(np int) {
	doParallel(func(np int) {
		// CalculateEdgeEulerFlux is where the Riemann problem is solved at the
		// neighbor faces, and the edge boundary conditions are applied.
		c.CalculateEdgeEulerFlux(rk.Time, rk.Q_Face, rk.Flux_Face,
			rk.EdgeQ1[np], rk.EdgeQ2[np], c.SortedEdgeKeys[np]) // Global
	})
	doParallel(func(np int) {
		if sd != nil {
			sd.CalculateEpsilonGradient(c, np, QQQAll[np])
		}
	})
	// doSerial(func(np int) {
	doParallel(func(np int) {
		c.StoreEdgeAggregates(rk.Epsilon, rk.Jdet, rk.Q_Face, c.SortedEdgeKeys[np])
		if sd != nil {
			c.StoreEdgeViscousFlux(rk.Epsilon, rk.EdgeQ1[np], c.SortedEdgeKeys[np])
		}
	})
	// doSerial(func(np int) {
	doParallel(func(np int) {
		if rkStep == 0 {
			c.InitializeDT(rk.DT[np])
		}
		rk.GlobalMaxWaveSpeed[np], rk.GlobalMaxViscWaveSpeed[np] =
			c.CalcElementMaxWaveSpeed(np, rk.DT[np], rk.DTVisc[np])
	})
	if !c.LocalTimeStepping {
		rk.calculateGlobalDT(c) // Compute the global DT for nonlocal timestepping - must be done serially
	}
	doParallel(func(np int) {
		c.CalculateLocalDT(rk.DT[np], rk.DTVisc[np], rk.GlobalDT)
		// Perform a Runge Kutta pseudo time step
		rkAdvance(np)
	})
	return
}

func (c *Euler) InitializeDT(DT utils.Vector) {
	var (
		KMax, _ = DT.Dims()
	)
	// Initialize local time stepping
	for k := 0; k < KMax; k++ {
		DT.DataP[k] = -100 // Global
	}
}

func (c *Euler) RHSInternalPoints(Kmax int, Jdet utils.Matrix, F_RT_DOF, RHSQ [4]utils.Matrix) {
	var (
		JdetD = Jdet.DataP
		Nint  = c.DFR.FluxElement.NpInt
		data  []float64
	)
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
	for n := 0; n < 4; n++ { // For each of the 4 equations in [rho, rhoU, rhoV, E]
		// Unit triangle divergence matrix times the flux projected onto the RT elements: F_RT_DOF => Div(Flux) in (r,s)
		c.DFR.FluxElement.DivInt.Mul(F_RT_DOF[n], RHSQ[n])
		// Multiply each element's divergence by 1/||J|| to go (r,s)->(x,y), and -1 for the RHS
		data = RHSQ[n].DataP
		for k := 0; k < Kmax; k++ {
			var (
				oojd = 1. / JdetD[k]
			)
			for i := 0; i < Nint; i++ {
				ind := k + i*Kmax
				// data[ind] /= -JdetD[k]
				data[ind] *= -oojd
			}
		}
	}
}

func (c *Euler) SetRTFluxInternal(Kmax int, Jdet, Jinv utils.Matrix, F_RT_DOF, Q [4]utils.Matrix) {
	var (
		Nint   = c.DFR.FluxElement.NpInt
		NpFlux = c.DFR.FluxElement.Np
		fdofD  = [4][]float64{F_RT_DOF[0].DataP, F_RT_DOF[1].DataP, F_RT_DOF[2].DataP, F_RT_DOF[3].DataP}
	)
	/*
		Zero out DOF storage to promote easier bug avoidance
	*/
	for n := 0; n < 4; n++ {
		for i := 0; i < Kmax*NpFlux; i++ {
			fdofD[n][i] = 0.
		}
	}
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
	case SHOCKTUBE:
		c.SolutionX = c.ShardByK(c.DFR.SolutionX)
		c.SolutionY = c.ShardByK(c.DFR.SolutionY)
		gamma := 1.4
		c.FSIn = NewFreestreamFromQinf(gamma, [4]float64{1, 0, 0, 1 / (gamma - 1)})
		c.FSOut = NewFreestreamFromQinf(gamma, [4]float64{0.125, 0, 0, 0.1 / (gamma - 1)})
		NP := c.Partitions.ParallelDegree
		for np := 0; np < NP; np++ {
			var (
				Kmax = c.Partitions.GetBucketDimension(np)
				Nint = c.DFR.SolutionElement.Np
			)
			for n := 0; n < 4; n++ {
				c.Q[np][n] = utils.NewMatrix(Nint, Kmax)
			}
			for k := 0; k < Kmax; k++ {
				for i := 0; i < Nint; i++ {
					ind := k + i*Kmax
					if c.SolutionX[np].DataP[ind] < 0.5 { // The length of the domain should be 0->1
						c.Q[np][0].DataP[ind] = c.FSIn.Qinf[0]
						c.Q[np][1].DataP[ind] = c.FSIn.Qinf[1]
						c.Q[np][2].DataP[ind] = c.FSIn.Qinf[2]
						c.Q[np][3].DataP[ind] = c.FSIn.Qinf[3]
					} else {
						c.Q[np][0].DataP[ind] = c.FSOut.Qinf[0]
						c.Q[np][1].DataP[ind] = c.FSOut.Qinf[1]
						c.Q[np][2].DataP[ind] = c.FSOut.Qinf[2]
						c.Q[np][3].DataP[ind] = c.FSOut.Qinf[3]
					}
				}
			}
		}
		// There are 5 points vertically, then we double sample (4 pts per cell) to capture all the dynamics
		c.ShockTube = sod_shock_tube.NewSODShockTube(4*c.DFR.K/5, c.DFR)
	case IVORTEX:
		c.FSFar = NewFreestreamFromQinf(1.4, [4]float64{1, 1, 0, 3})
		c.SolutionX = c.ShardByK(c.DFR.SolutionX)
		c.SolutionY = c.ShardByK(c.DFR.SolutionY)
		NP := c.Partitions.ParallelDegree
		for np := 0; np < NP; np++ {
			c.AnalyticSolution, c.Q[np] = c.InitializeIVortex(c.SolutionX[np], c.SolutionY[np])
		}
		// SetScalar "Wall" BCs to IVortex
		var count int
		for _, e := range c.DFR.Tris.Edges {
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

func (c *Euler) CheckIfFinished(Time float64, steps int) (finished bool) {
	if Time >= c.FinalTime || steps >= c.MaxIterations {
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
		fmt.Printf("    iter    time      dt")
	}
	fmt.Printf("       Res0       Res1       Res2")
	fmt.Printf("       Res3         L1         L2\n")
}
func (c *Euler) PrintUpdate(Time, dt float64, steps int, Q, Residual [][4]utils.Matrix, printMem bool, limitedPoints []int) {
	format := "%11.4e"
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
	var lpsum int
	for _, val := range limitedPoints {
		lpsum += val
	}
	if lpsum != 0 {
		fmt.Printf(" #limited:%5d/%-8d ", lpsum, c.DFR.K*c.DFR.SolutionElement.Np)
	}
	fmt.Printf("\n")
}

func (c *Euler) PrintFinal(elapsed time.Duration, steps int) {
	rate := float64(elapsed.Microseconds()) / (float64(c.DFR.K * steps))
	fmt.Printf("\nRate of execution = %8.5f us/(element*iteration) over %d iterations\n", rate, steps)
	instructionRate := 1344.282 // measured instructions per element per iteration
	iRate := instructionRate * float64(c.DFR.K*steps) / elapsed.Seconds()
	var mem syscall.Rusage
	_ = syscall.Getrusage(syscall.RUSAGE_SELF, &mem)
	usec, ssec := float64(mem.Utime.Nano())/1.e9, float64(mem.Stime.Nano())/1.e9
	fmt.Printf("OS Stats: Utime=%6.2f(s), Stime=%6.2f(s), RSS Memory=%6.2f(MB)\n", usec, ssec, float64(mem.Maxrss)/1024.)
	fmt.Printf("Estimated op rate = %8.5f Giga-ops / second, based on measure op rate of %8.5f per element per iteration\n",
		iRate/1000000000., instructionRate)
}

func (c *Euler) GetSolutionGradientUsingRTElement(myThread, varNum int, Q [4]utils.Matrix, GradX, GradY, DOFX, DOFY utils.Matrix) {
	/*
		Dimensions:
			Q[4] should be NpInt x Kmax
			All others should be NpFlux x Kmax
	*/
	var (
		dfr          = c.DFR
		NpInt        = dfr.FluxElement.NpInt
		NpEdge       = dfr.FluxElement.NpEdge
		KmaxGlobal   = c.Partitions.MaxIndex
		Kmax         = c.Partitions.GetBucketDimension(myThread)
		DXmd, DYmd   = dfr.DXMetric.DataP, dfr.DYMetric.DataP
		DOFXd, DOFYd = DOFX.DataP, DOFY.DataP
	)
	var Un float64
	for k := 0; k < Kmax; k++ {
		kGlobal := c.Partitions.GetGlobalK(k, myThread)
		for i := 0; i < NpInt; i++ {
			ind := k + i*Kmax
			indG := kGlobal + i*KmaxGlobal
			ind2 := k + (i+NpInt)*Kmax
			ind2G := kGlobal + (i+NpInt)*KmaxGlobal
			Un = Q[varNum].DataP[ind]
			DOFXd[ind], DOFYd[ind] = DXmd[indG]*Un, DYmd[indG]*Un
			DOFXd[ind2], DOFYd[ind2] = DXmd[ind2G]*Un, DYmd[ind2G]*Un
		}
	}
	// Load edge solution values from solution storage
	for k := 0; k < Kmax; k++ {
		kGlobal := c.Partitions.GetGlobalK(k, myThread)
		for edgeNum := 0; edgeNum < 3; edgeNum++ {
			edgeVals, sign := c.EdgeStore.GetEdgeValues(EdgeQValues, myThread, k, varNum, edgeNum, dfr)
			var (
				ii    int
				shift = edgeNum * NpEdge
			)
			for i := 0; i < NpEdge; i++ {
				if sign < 0 {
					ii = NpEdge - 1 - i
				} else {
					ii = i
				}
				ind := k + (2*NpInt+i+shift)*Kmax
				indG := kGlobal + (2*NpInt+i+shift)*KmaxGlobal
				Un = edgeVals[ii]
				DOFXd[ind] = DXmd[indG] * Un
				DOFYd[ind] = DYmd[indG] * Un
			}
		}
	}
	// Calculate Grad(U)
	dfr.FluxElement.Div.Mul(DOFX, GradX) // R Derivative, divergence x RT_DOF is R derivative for this DOF
	dfr.FluxElement.Div.Mul(DOFY, GradY) // S Derivative, divergence x RT_DOF is S derivative for this DOF
}

func (c *Euler) GetSolutionGradient(myThread, varNum int, Q [4]utils.Matrix, GradX, GradY, DR, DS utils.Matrix) {
	/*
		Dimensions:
			Q[4] should be NpInt x Kmax
			All others should be NpFlux x Kmax
	*/
	var (
		dfr    = c.DFR
		NpFlux = dfr.FluxElement.Np
		Kmax   = c.Partitions.GetBucketDimension(myThread)
	)
	dfr.FluxDr.Mul(Q[varNum], DR)
	dfr.FluxDs.Mul(Q[varNum], DS)
	for k := 0; k < Kmax; k++ {
		var (
			JinvD = dfr.Jinv.DataP[4*k : 4*(k+1)]
		)
		for i := 0; i < NpFlux; i++ {
			ind := k + i*Kmax
			GradX.DataP[ind] = JinvD[0]*DR.DataP[ind] + JinvD[2]*DS.DataP[ind]
			GradY.DataP[ind] = JinvD[1]*DR.DataP[ind] + JinvD[3]*DS.DataP[ind]
		}
	}
}

func (rk *RungeKutta5SSP) calculateGlobalDT(c *Euler) {
	var (
		pm = c.Partitions
		NP = pm.ParallelDegree
	)
	// Must be done in controller/sync process
	var globalMaxWaveSpeed float64
	for np := 0; np < NP; np++ {
		globalMaxWaveSpeed = max(globalMaxWaveSpeed, rk.GlobalMaxWaveSpeed[np])
	}
	rk.GlobalDT = c.CFL / globalMaxWaveSpeed
	if c.Dissipation != nil {
		var (
			NP1    = float64(c.DFR.N + 1)
			NP12   = NP1 * NP1
			C_diff = 1. / NP12
		)
		var globalMaxViscWaveSpeed float64
		for np := 0; np < NP; np++ {
			globalMaxViscWaveSpeed = max(globalMaxViscWaveSpeed, rk.GlobalMaxViscWaveSpeed[np])
		}
		rk.GlobalDT = min(rk.GlobalDT, C_diff/globalMaxViscWaveSpeed)
	}
	if rk.Time+rk.GlobalDT > c.FinalTime {
		rk.GlobalDT = c.FinalTime - rk.Time
	}
}

func (c *Euler) CalculateLocalDT(DT, DTVisc utils.Vector, GlobalDT float64) {
	var (
		Kmax, _ = DT.Dims()
		NP1     = float64(c.DFR.N + 1)
		NP12    = NP1 * NP1
	)
	// c.LocalTimeStepping = false
	if !c.LocalTimeStepping {
		for k := range DT.DataP {
			DT.DataP[k] = GlobalDT
		}
	} else {
		for k := 0; k < Kmax; k++ {
			// Set each element's DT to CFL*h/(max_wave_speed)
			DT.DataP[k] = c.CFL / DT.DataP[k]
		}
		if c.Dissipation != nil {
			AVThreshold := 1.e-9
			// C_diff≈0.1–0.25;
			// C_diff := 0.25
			C_diff := 1. / NP12
			for k := 0; k < Kmax; k++ {
				if DTVisc.DataP[k] > AVThreshold {
					DTVisc.DataP[k] = C_diff / DTVisc.DataP[k]
					DT.DataP[k] = min(DT.DataP[k], DTVisc.DataP[k])
				}
			}
		}
	}
}

func (c *Euler) UpdateElementMean(Q [4]utils.Matrix, QMean [4]utils.Vector) {
	var (
		Np      = c.DFR.SolutionElement.Np
		_, KMax = Q[0].Dims()
		P       = c.DFR.N
	)
	for n := 0; n < 4; n++ {
		for k := 0; k < KMax; k++ {
			QMean[n].DataP[k] = 0.
			for i := 0; i < Np; i++ {
				QMean[n].DataP[k] +=
					Q[n].At(i, k) * DG2D.WilliamsShunnJamesonCoords[P][i].W
			}
		}
	}
	return
}

func (c *Euler) shardEtoV(EtoV utils.Matrix) (ev []utils.Matrix) {
	var (
		pm = c.Partitions
		NP = pm.ParallelDegree
		// KMax, _ = EtoV.Dims()
	)
	ev = make([]utils.Matrix, NP)
	for np := 0; np < NP; np++ {
		KMax := pm.GetBucketDimension(np)
		ev[np] = utils.NewMatrix(KMax, 3)
		kmin, kmax := pm.GetBucketRange(np)
		var klocal int
		for k := kmin; k < kmax; k++ {
			ev[np].Set(klocal, 0, EtoV.At(k, 0))
			ev[np].Set(klocal, 1, EtoV.At(k, 1))
			ev[np].Set(klocal, 2, EtoV.At(k, 2))
			klocal++
		}
		if klocal != KMax {
			msg := fmt.Errorf("dimension incorrect, should be %d, is %d", KMax, klocal)
			panic(msg)
		}
	}
	return
}

func (c *Euler) MergeElementScalarToVertices(VtoE VertexToElement,
	ElementScalar []utils.Vector, ElementVertexScalar []float64,
	mergeFunc func(a, b float64) float64) {
	// ElementScalar is the sharded scalar field defined on K elements
	// ElementVertexScalar is the target for consolidation of element values
	// mergeFunc is typically one of min/max/sum
	oldVert := -1
	for _, val := range VtoE {
		vert, k, threadID := int(val[0]), int(val[1]), int(val[2])
		if oldVert == vert { // we're in the middle of processing this vert, update normally
			ElementVertexScalar[vert] = mergeFunc(ElementVertexScalar[vert],
				ElementScalar[threadID].AtVec(k))
		} else { // we're on a new vertex, reset the vertex value
			ElementVertexScalar[vert] = ElementScalar[threadID].AtVec(k)
			oldVert = vert
		}
	}
}

func (c *Euler) MergeVertexScalarToElement(ElementVertexScalar []float64,
	ElementScalar utils.Vector,
	EtoV utils.Matrix,
	mergeFunc func(a, b float64) float64, finalFunc func(a float64) float64) {
	// ElementVertexScalar is the scalar value on each vertex, 3 per element
	// ElementScalar is the non-sharded target element
	var (
		KMax, _ = EtoV.Dims()
	)
	for k := 0; k < KMax; k++ {
		scalarAccum := 0.
		for v := 0; v < 3; v++ {
			ind := v + 3*k
			scalarAccum =
				mergeFunc(scalarAccum, ElementVertexScalar[int(EtoV.DataP[ind])])
		}
		ElementScalar.Set(k, finalFunc(scalarAccum))
	}
	return
}

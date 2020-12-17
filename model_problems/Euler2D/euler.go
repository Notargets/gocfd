package Euler2D

import (
	"fmt"
	"image/color"
	"math"
	"runtime"
	"sort"
	"strings"
	"sync"
	"time"

	"github.com/notargets/gocfd/types"

	graphics2D "github.com/notargets/avs/geometry"

	"github.com/notargets/avs/functions"

	"github.com/notargets/gocfd/model_problems/Euler2D/isentropic_vortex"

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
	Qinf                  [4]float64
	dfr                   *DG2D.DFR2D
	Q                     [4]utils.Matrix // Solution variables, stored at solution point locations, Np_solution x K
	DT                    utils.Matrix    // Local time step for steady state solution
	Q_Face                [4]utils.Matrix // Solution variables, interpolated to and stored at edge point locations, Np_edge x K
	F_RT_DOF              [4]utils.Matrix // Normal Projected Flux, stored at flux/solution point locations, Np_flux x K
	chart                 ChartState
	FluxCalcAlgo          FluxType
	Case                  InitType
	AnalyticSolution      ExactState
	FluxCalcMock          func(Gamma, rho, rhoU, rhoV, E float64) (Fx, Fy [4]float64) // For testing
	SortedEdgeKeys        EdgeKeySlice
	ParallelDegree        int // Number of go routines to use for parallel execution
	LocalTimeStepping     bool
}

type ChartState struct {
	chart *chart2d.Chart2D
	fs    *functions.FSurface
	gm    *graphics2D.TriMesh
}

type ExactState interface {
	GetStateC(t, x, y float64) (rho, rhoU, rhoV, E float64)
	GetDivergence(t, x, y float64) (div [4]float64)
}

type InitType uint

const (
	FREESTREAM InitType = iota
	IVORTEX
)

var (
	InitNames = map[string]InitType{
		"freestream": FREESTREAM,
		"ivortex":    IVORTEX,
	}
	InitPrintNames = []string{"Freestream", "Inviscid Vortex Analytic Solution"}
)

func NewInitType(label string) (it InitType) {
	var (
		ok  bool
		err error
	)
	if len(label) == 0 {
		err = fmt.Errorf("empty init type, must be one of %v", InitNames)
		panic(err)
	}
	label = strings.ToLower(label)
	if it, ok = InitNames[label]; !ok {
		err = fmt.Errorf("unable to use init type named %s", label)
		panic(err)
	}
	return
}

func (it InitType) Print() (txt string) {
	txt = InitPrintNames[it]
	return
}

type FluxType uint

const (
	FLUX_Average FluxType = iota
	FLUX_LaxFriedrichs
	FLUX_Roe
)

var (
	FluxNames = map[string]FluxType{
		"average": FLUX_Average,
		"lax":     FLUX_LaxFriedrichs,
		"roe":     FLUX_Roe,
	}
	FluxPrintNames = []string{"Average", "Lax Friedrichs", "Roe"}
)

func (ft FluxType) Print() (txt string) {
	txt = FluxPrintNames[ft]
	return
}

func NewFluxType(label string) (ft FluxType) {
	var (
		ok  bool
		err error
	)
	label = strings.ToLower(label)
	if ft, ok = FluxNames[label]; !ok {
		err = fmt.Errorf("unable to use flux named %s", label)
		panic(err)
	}
	return
}

func (c *Euler) CalcFS(Minf, Gamma, Alpha float64) {
	var (
		ooggm1 = 1. / (Gamma * (Gamma - 1.))
	)
	c.Qinf[0] = 1
	c.Qinf[1] = Minf * math.Cos(Alpha*math.Pi/180.)
	c.Qinf[2] = Minf * math.Sin(Alpha*math.Pi/180.)
	c.Qinf[3] = ooggm1 + 0.5*Minf*Minf
}

func NewEuler(FinalTime float64, N int, meshFile string, CFL float64, fluxType FluxType, Case InitType, ProcLimit int, Minf, Gamma, Alpha float64, LocalTime, plotMesh, verbose bool) (c *Euler) {
	c = &Euler{
		MeshFile:          meshFile,
		CFL:               CFL,
		FinalTime:         FinalTime,
		FluxCalcAlgo:      fluxType,
		Case:              Case,
		Gamma:             Gamma,
		FluxCalcMock:      FluxCalc,
		LocalTimeStepping: LocalTime,
	}
	if ProcLimit != 0 {
		c.ParallelDegree = ProcLimit
	} else {
		c.ParallelDegree = runtime.NumCPU()
	}
	runtime.GOMAXPROCS(runtime.NumCPU())

	c.dfr = DG2D.NewDFR2D(N, plotMesh, meshFile)
	c.InitializeMemory()
	if verbose {
		fmt.Printf("Euler Equations in 2 Dimensions\n")
		fmt.Printf("Using %d go routines in parallel\n", runtime.NumCPU())
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
	/*
		fmt.Printf("Total number of edges = %d\n", 3*c.dfr.K)
		var internalEdgesCount, boundaryEdgesCount, unconnectedEdgesCount int
		for _, e := range c.dfr.Tris.Edges {
			switch e.NumConnectedTris {
			case 0:
				unconnectedEdgesCount++
			case 1:
				boundaryEdgesCount++
			case 2:
				internalEdgesCount++
			}
		}
		fmt.Printf("Count of internal edges, incl Periodic Boundaries = %d\n", internalEdgesCount)
		fmt.Printf("Count of unhandled edges = %d\n", unconnectedEdgesCount)
		fmt.Printf("Count of boundary edges = %d\n", boundaryEdgesCount)
	*/
	switch c.Case {
	case FREESTREAM:
		c.CalcFS(Minf, Gamma, Alpha)
		c.InitializeFS()
	case IVORTEX:
		c.Qinf = [4]float64{1, 1, 0, 3}
		c.AnalyticSolution, c.Q = c.InitializeIVortex(c.dfr.SolutionX, c.dfr.SolutionY)
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
		fmt.Printf("done\n")
		fmt.Printf("Algorithm: %s\n", c.FluxCalcAlgo.Print())
		fmt.Printf("CFL = %8.4f, Polynomial Degree N = %d (1 is linear), Num Elements K = %d\n\n\n", CFL, N, c.dfr.K)
	}
	return
}

func (c *Euler) Solve(pm *PlotMeta) {
	var (
		FinalTime        = c.FinalTime
		Time, dt         float64
		Q1, Q2, Residual [4]utils.Matrix
		steps            int
		Np               = c.dfr.SolutionElement.Np
		Kmax             = c.dfr.K
		plotQ            = pm.Plot
	)
	// Initialize memory for RHS
	for n := 0; n < 4; n++ {
		Q1[n] = utils.NewMatrix(Np, Kmax)
		Q2[n] = utils.NewMatrix(Np, Kmax)
		Residual[n] = Q1[n] // optimize memory using an alias
	}
	fmt.Printf("solving until finaltime = %8.5f\n", FinalTime)
	c.InterpolateSolutionToEdges(c.Q)
	start := time.Now()
	for {
		if Time >= FinalTime {
			break
		}
		dt = c.CalculateDT()
		if Time+dt > FinalTime {
			dt = FinalTime - Time
		}
		rhsQ := c.RHS(c.Q, Time)
		for n := 0; n < 4; n++ {
			if c.LocalTimeStepping {
				Q1[n].Apply4Parallel(rhsQ[n], c.Q[n], c.DT, func(q1, rhsq, q, dT float64) (res float64) {
					res = q + rhsq*dT
					return
				}, c.ParallelDegree)
			} else {
				Q1[n].Apply3Parallel(rhsQ[n], c.Q[n], func(q1, rhsq, q float64) (res float64) {
					res = q + rhsq*dt
					return
				}, c.ParallelDegree)
			}
		}
		rhsQ = c.RHS(Q1, Time)
		for n := 0; n < 4; n++ {
			if c.LocalTimeStepping {
				Q2[n].Apply5Parallel(rhsQ[n], Q1[n], c.Q[n], c.DT, func(q2, rhsq, q1, q, dT float64) (res float64) {
					res = 0.25 * (q1 + 3*q + rhsq*dT)
					return
				}, c.ParallelDegree)
			} else {
				Q2[n].Apply4Parallel(rhsQ[n], Q1[n], c.Q[n], func(q2, rhsq, q1, q float64) (res float64) {
					res = 0.25 * (q1 + 3*q + rhsq*dt)
					return
				}, c.ParallelDegree)
			}
		}
		rhsQ = c.RHS(Q2, Time)
		for n := 0; n < 4; n++ {
			if c.LocalTimeStepping {
				Residual[n].Apply5Parallel(rhsQ[n], Q2[n], c.Q[n], c.DT, func(resid, rhsq, q2, q, dT float64) (res float64) {
					res = (2. / 3) * (q2 - q + dT*rhsq)
					return
				}, c.ParallelDegree)
			} else {
				Residual[n].Apply4Parallel(rhsQ[n], Q2[n], c.Q[n], func(resid, rhsq, q2, q float64) (res float64) {
					res = (2. / 3) * (q2 - q + dt*rhsq)
					return
				}, c.ParallelDegree)
			}
			c.Q[n].Add(Residual[n])
		}
		steps++
		Time += dt
		if plotQ && steps%pm.StepsBeforePlot == 0 {
			c.PlotQ(c.Q, pm) // wait till we implement time iterative frame updates
		}
		if steps%pm.StepsBeforePlot == 0 {
			fmt.Printf("\nTime,dt = %8.5f,%8.5f, Residual[eq#]Min/Max:", Time, dt)
			/*
				var Mmin, Mmax float64
				Mmax, Mmin = -1, 1000
				for ik := 0; ik < c.dfr.K*c.dfr.SolutionElement.Np; ik++ {
					C, u, v := GetSpeedsDirect(ik, c.Gamma, c.Q)
					U := math.Sqrt(u*u + v*v)
					M := U / C
					if M < Mmin {
						Mmin = M
					}
					if M > Mmax {
						Mmax = M
					}
				}
				fmt.Printf(" Mach min:%8.5f max:%8.5f ", Mmin, Mmax)
			*/
			for n := 0; n < 4; n++ {
				fmt.Printf(" [%d] %8.5f,%8.5f ", n, Residual[n].Min(), Residual[n].Max())
			}
		}
	}
	elapsed := time.Now().Sub(start)
	rate := float64(elapsed.Microseconds()) / (float64(Kmax * steps))
	fmt.Printf("\nRate of execution = %8.5f us/(element*iteration) over %d iterations\n", rate, steps)
}

func (c *Euler) PlotQ(Q [4]utils.Matrix, pm *PlotMeta) {
	var (
		plotField = pm.Field
		delay     = pm.FrameTime
		lineType  = pm.LineType
		scale     = pm.Scale
		translate = [2]float32{float32(pm.TranslateX), float32(pm.TranslateY)}
		oField    = c.GetPlotField(Q, plotField)
		fI        = c.dfr.ConvertScalarToOutputMesh(oField)
	)

	if c.chart.gm == nil {
		c.chart.gm = c.dfr.OutputMesh()
	}
	c.chart.fs = functions.NewFSurface(c.chart.gm, [][]float32{fI}, 0)
	fmt.Printf(" Plot>%s min,max = %8.5f,%8.5f\n", pm.Field.String(), oField.Min(), oField.Max())
	c.PlotFS(pm.FieldMinP, pm.FieldMaxP, 0.95*oField.Min(), 1.05*oField.Max(), scale, translate, lineType)
	utils.SleepFor(int(delay.Milliseconds()))
	return
}

func (c *Euler) GetPlotField(Q [4]utils.Matrix, plotField PlotField) (field utils.Matrix) {
	var (
		Kmax = c.dfr.K
		Np   = c.dfr.SolutionElement.Np
	)
	switch plotField {
	case Density, XMomentum, YMomentum, Energy:
		field = c.dfr.FluxInterpMatrix.Mul(Q[int(plotField)])
	case Mach:
		fld := utils.NewMatrix(Np, Kmax)
		for ik := 0; ik < Kmax*Np; ik++ {
			C, u, v := GetSpeedsDirect(ik, c.Gamma, Q)
			fld.Data()[ik] = math.Sqrt(u*u+v*v) / C
		}
		field = c.dfr.FluxInterpMatrix.Mul(fld)
	}
	return
}

func (c *Euler) PlotFS(fminP, fmaxP *float64, fmin, fmax float64, scale float64, translate [2]float32, ltO ...chart2d.LineType) {
	var (
		fs      = c.chart.fs
		trimesh = fs.Tris
		lt      = chart2d.NoLine
	)
	if c.chart.chart == nil {
		box := graphics2D.NewBoundingBox(trimesh.GetGeometry())
		box = box.Scale(float32(scale))
		box = box.Translate(translate)
		c.chart.chart = chart2d.NewChart2D(1900, 1080, box.XMin[0], box.XMax[0], box.XMin[1], box.XMax[1])
		colorMap := utils2.NewColorMap(float32(fmin), float32(fmax), 1.)
		c.chart.chart.AddColorMap(colorMap)
		go c.chart.chart.Plot()
	}

	if fminP != nil || fmaxP != nil {
		switch {
		case fminP != nil && fmaxP != nil:
			fmin, fmax = *fminP, *fmaxP
		case fminP != nil:
			fmin = *fminP
		case fmaxP != nil:
			fmax = *fmaxP
		}
		colorMap := utils2.NewColorMap(float32(fmin), float32(fmax), 1.)
		c.chart.chart.AddColorMap(colorMap)
	}

	white := color.RGBA{R: 255, G: 255, B: 255, A: 1}
	black := color.RGBA{R: 0, G: 0, B: 0, A: 1}
	_, _ = white, black
	if len(ltO) != 0 {
		lt = ltO[0]
	}
	if err := c.chart.chart.AddFunctionSurface("FSurface", *fs, lt, white); err != nil {
		panic("unable to add function surface series")
	}
}

func (c *Euler) CalculateDT() (dt float64) {
	var (
		Np1      = c.dfr.N + 1
		Np12     = float64(Np1 * Np1)
		wsMaxAll = -math.MaxFloat64
		wsMax    = make([]float64, c.ParallelDegree)
		wg       = sync.WaitGroup{}
	)
	// Setup max wavespeed before loop
	for k := 0; k < c.dfr.K; k++ {
		c.DT.Data()[k] = -100
	}
	for nn := 0; nn < c.ParallelDegree; nn++ {
		wsMax[nn] = wsMaxAll
	}
	// Loop over all edges, calculating max wavespeed
	for nn := 0; nn < c.ParallelDegree; nn++ {
		ind, end := c.split1D(len(c.SortedEdgeKeys), nn)
		wg.Add(1)
		go func(ind, end, nn int) {
			for ii := ind; ii < end; ii++ {
				edgeKey := c.SortedEdgeKeys[ii]
				e := c.dfr.Tris.Edges[edgeKey]
				// for _, e := range c.dfr.Tris.Edges {
				var (
					edgeLen = e.GetEdgeLength()
					Nedge   = c.dfr.FluxElement.Nedge
				)
				conn := 0
				var (
					k       = int(e.ConnectedTris[conn])
					edgeNum = int(e.ConnectedTriEdgeNumber[conn])
					shift   = edgeNum * Nedge
				)
				Jdet := c.dfr.Jdet.At(k, 0)
				// fmt.Printf("N, Np12, edgelen, Jdet = %d,%8.5f,%8.5f,%8.5f\n", c.dfr.N, Np12, edgeLen, Jdet)
				fs := 0.5 * Np12 * edgeLen / Jdet
				edgeMax := -100.
				for i := shift; i < shift+Nedge; i++ {
					C, u, v := c.GetSpeeds(k, i, c.Q_Face)
					waveSpeed := fs * (math.Sqrt(u*u+v*v) + C)
					wsMax[nn] = math.Max(waveSpeed, wsMax[nn])
					if waveSpeed > edgeMax {
						edgeMax = waveSpeed
					}
				}
				if edgeMax > c.DT.Data()[k] {
					c.DT.Data()[k] = edgeMax
				}
				if e.NumConnectedTris == 2 { // Add the wavespeed to the other tri connected to this edge if needed
					k = int(e.ConnectedTris[1])
					if edgeMax > c.DT.Data()[k] {
						c.DT.Data()[k] = edgeMax
					}
				}
			}
			wg.Done()
		}(ind, end, nn)
	}
	wg.Wait()
	// Replicate local time step to the other solution points for each k
	for k := 0; k < c.dfr.K; k++ {
		c.DT.Data()[k] = c.CFL / c.DT.Data()[k]
	}
	for i := 1; i < c.dfr.SolutionElement.Np; i++ {
		for k := 0; k < c.dfr.K; k++ {
			ind := k + c.dfr.K*i
			c.DT.Data()[ind] = c.DT.Data()[k]
		}
	}
	for nn := 0; nn < c.ParallelDegree; nn++ {
		wsMaxAll = math.Max(wsMaxAll, wsMax[nn])
	}
	dt = c.CFL / wsMaxAll
	return
}

func (c *Euler) GetState(k, i int, Q [4]utils.Matrix) (rho, u, v, p, C, E float64) {
	var (
		Gamma = c.Gamma
		Kmax  = c.dfr.K
		ind   = k + Kmax*i
		q     = [4]float64{Q[0].Data()[ind], Q[1].Data()[ind], Q[2].Data()[ind], Q[3].Data()[ind]}
	)
	rho, u, v, p, C, E = GetPrimitiveVariables(Gamma, q)
	return
}

func (c *Euler) GetPressure(k, i int, Q [4]utils.Matrix) (p float64) {
	var (
		Gamma              = c.Gamma
		GM1                = Gamma - 1.
		Kmax               = c.dfr.K
		ind                = k + Kmax*i
		rho, rhou, rhov, E = Q[0].Data()[ind], Q[1].Data()[ind], Q[2].Data()[ind], Q[3].Data()[ind]
	)
	p = GM1 * (E - 0.5*(rhou*rhou+rhov*rhov)/rho)
	return
}

func (c *Euler) GetSpeeds(k, i int, Q [4]utils.Matrix) (C, u, v float64) {
	var (
		Kmax = c.dfr.K
		ind  = k + Kmax*i
	)
	C, u, v = GetSpeedsDirect(ind, c.Gamma, Q)
	return
}

func GetSpeedsDirect(ind int, gamma float64, Q [4]utils.Matrix) (C, u, v float64) {
	var (
		Gamma              = gamma
		GM1                = Gamma - 1.
		rho, rhou, rhov, E = Q[0].Data()[ind], Q[1].Data()[ind], Q[2].Data()[ind], Q[3].Data()[ind]
		oorho              = 1. / rho
		p                  = GM1 * (E - 0.5*(rhou*rhou+rhov*rhov)*oorho)
	)
	C = math.Sqrt(math.Abs(Gamma * p * oorho))
	u = rhou * oorho
	v = rhov * oorho
	return
}

func GetPrimitiveVariables(Gamma float64, Q [4]float64) (rho, u, v, p, C, E float64) {
	var (
		GM1 = Gamma - 1.
	)
	rho = Q[0]
	u = Q[1] / rho
	v = Q[2] / rho
	E = Q[3]
	p = GM1 * (E - 0.5*rho*(u*u+v*v))
	C = math.Sqrt(math.Abs(Gamma * p / rho))
	return
}

func (c *Euler) RHS(Q [4]utils.Matrix, Time float64) (RHSCalc [4]utils.Matrix) {
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
	c.AssembleRTNormalFlux(Q, Time) // Assembles F_RT_DOF for use in calculations using RT element
	var wg = sync.WaitGroup{}
	for n := 0; n < 4; n++ {
		wg.Add(1)
		go func(n int) {
			RHSCalc[n] = c.dfr.FluxElement.DivInt.Mul(c.F_RT_DOF[n]) // Calculate divergence for the internal node points
			c.DivideByJacobian(c.dfr.FluxElement.Nint, RHSCalc[n].Data(), -1)
			wg.Done()
		}(n)
	}
	wg.Wait()
	return
}

func (c *Euler) DivideByJacobian(Nmax int, data []float64, scale float64) {
	var (
		Kmax = c.dfr.K
	)
	for k := 0; k < Kmax; k++ {
		Jdet := c.dfr.Jdet.At(k, 0)
		for i := 0; i < Nmax; i++ {
			ind := k + i*Kmax
			data[ind] /= (Jdet * scale) // Multiply divergence by -1 to produce the RHS
		}
	}
}

func (c *Euler) AssembleRTNormalFlux(Q [4]utils.Matrix, Time float64) {
	/*
		Solver approach:
		0) Solution is stored on sol points as Q
		0a) Flux is computed and stored in X, Y component projections in the 2*Nint front of F_RT_DOF
		1) Solution is extrapolated to edge points in Q_Face from Q
		2) Edges are traversed, flux is calculated and projected onto edge face normals, scaled and placed into F_RT_DOF
	*/
	/*
		Zero out DOF storage and Q_Face to promote easier bug avoidance
	*/
	for n := 0; n < 4; n++ {
		c.Q_Face[n].Scale(0.)
		c.F_RT_DOF[n].Scale(0.)
	}
	c.SetNormalFluxInternal(Q)           // Updates F_RT_DOF with values from Q
	c.InterpolateSolutionToEdges(Q)      // Interpolates Q_Face values from Q
	c.ParallelSetNormalFluxOnEdges(Time) // Updates F_RT_DOG with values from edges, including BCs and connected tris
}

func (c *Euler) SetNormalFluxInternal(Q [4]utils.Matrix) {
	var (
		Kmax = c.dfr.K
		Nint = c.dfr.FluxElement.Nint
		wg   = sync.WaitGroup{}
	)
	// Calculate flux and project into R and S (transformed) directions for the internal points
	//for k := 0; k < Kmax; k++ {
	for np := 0; np < c.ParallelDegree; np++ {
		ind, end := c.split1D(Kmax, np)
		wg.Add(1)
		go func() {
			for k := ind; k < end; k++ {
				for i := 0; i < Nint; i++ {
					ind := k + i*Kmax
					ind2 := k + (i+Nint)*Kmax
					Fr, Fs := c.CalculateFluxTransformed(k, i, Q)
					for n := 0; n < 4; n++ {
						rtD := c.F_RT_DOF[n].Data()
						rtD[ind], rtD[ind2] = Fr[n], Fs[n]
					}
				}
			}
			wg.Done()
		}()
	}
	wg.Wait()
}

func (c *Euler) InterpolateSolutionToEdges(Q [4]utils.Matrix) {
	// Interpolate from solution points to edges using precomputed interpolation matrix
	var wg = sync.WaitGroup{}
	for n := 0; n < 4; n++ {
		wg.Add(1)
		go func(n int) {
			c.Q_Face[n] = c.dfr.FluxEdgeInterpMatrix.Mul(Q[n])
			wg.Done()
		}(n)
	}
	wg.Wait()
}

type EdgeKeySlice []types.EdgeKey

func (p EdgeKeySlice) Len() int           { return len(p) }
func (p EdgeKeySlice) Less(i, j int) bool { return p[i] < p[j] }
func (p EdgeKeySlice) Swap(i, j int)      { p[i], p[j] = p[j], p[i] }

// Sort is a convenience method.
func (p EdgeKeySlice) Sort() { sort.Sort(p) }

func (c *Euler) ParallelSetNormalFluxOnEdges(Time float64) {
	var (
		Ntot = len(c.SortedEdgeKeys)
		wg   = sync.WaitGroup{}
	)
	var ind, end int
	for n := 0; n < c.ParallelDegree; n++ {
		ind, end = c.split1D(Ntot, n)
		wg.Add(1)
		go func(ind, end int) {
			c.SetNormalFluxOnEdges(Time, c.SortedEdgeKeys[ind:end])
			wg.Done()
		}(ind, end)
	}
	wg.Wait()
}

func (c *Euler) SetNormalFluxOnEdges(Time float64, edgeKeys EdgeKeySlice) {
	var (
		dfr                            = c.dfr
		Nint                           = dfr.FluxElement.Nint
		Nedge                          = dfr.FluxElement.Nedge
		Kmax                           = dfr.K
		normalFlux, normalFluxReversed = make([][4]float64, Nedge), make([][4]float64, Nedge)
	)
	for _, en := range edgeKeys {
		e := dfr.Tris.Edges[en]
		//for en, e := range dfr.Tris.Edges {
		switch e.NumConnectedTris {
		case 0:
			panic("unable to handle unconnected edges")
		case 1: // Handle edges with only one triangle - default is edge flux, which will be replaced by a BC flux
			var (
				k           = int(e.ConnectedTris[0])
				edgeNumber  = int(e.ConnectedTriEdgeNumber[0])
				shift       = edgeNumber * Nedge
				riemann     = true
				processFlux bool
			)
			processFlux = true
			normal, _ := c.getEdgeNormal(0, e, en)
			switch e.BCType {
			case types.BC_Far:
				for i := 0; i < Nedge; i++ {
					iL := i + shift
					ind := k + iL*Kmax
					QBC := c.RiemannBC(k, iL, c.Qinf, normal)
					c.Q_Face[0].Data()[ind] = QBC[0]
					c.Q_Face[1].Data()[ind] = QBC[1]
					c.Q_Face[2].Data()[ind] = QBC[2]
					c.Q_Face[3].Data()[ind] = QBC[3]
				}
			case types.BC_IVortex:
				// Set the flow variables to the exact solution
				X, Y := c.dfr.FluxX.Data(), c.dfr.FluxY.Data()
				for i := 0; i < Nedge; i++ {
					iL := i + shift
					indFlux := k + (2*Nint+iL)*Kmax
					x, y := X[indFlux], Y[indFlux]
					iv := c.AnalyticSolution.(*isentropic_vortex.IVortex)
					rho, rhoU, rhoV, E := iv.GetStateC(Time, x, y)
					ind := k + iL*Kmax
					var QBC [4]float64
					if riemann {
						Qinf := [4]float64{rho, rhoU, rhoV, E}
						QBC = c.RiemannBC(k, iL, Qinf, normal)
					} else {
						QBC = [4]float64{rho, rhoU, rhoV, E}
					}
					c.Q_Face[0].Data()[ind] = QBC[0]
					c.Q_Face[1].Data()[ind] = QBC[1]
					c.Q_Face[2].Data()[ind] = QBC[2]
					c.Q_Face[3].Data()[ind] = QBC[3]
				}
			case types.BC_Wall, types.BC_Cyl:
				processFlux = false
				for i := 0; i < Nedge; i++ {
					ie := i + shift
					p := c.GetPressure(k, ie, c.Q_Face) // Get pressure
					for n := 0; n < 4; n++ {
						normalFlux[i][0] = 0
						normalFlux[i][3] = 0
						normalFlux[i][1] = normal[0] * p
						normalFlux[i][2] = normal[1] * p
					}
				}
				c.SetNormalFluxOnRTEdge(k, edgeNumber, normalFlux, e.IInII[0])
			case types.BC_PeriodicReversed, types.BC_Periodic:
				processFlux = false
			}
			if processFlux {
				var Fx, Fy [4]float64
				for i := 0; i < Nedge; i++ {
					ie := i + shift
					Fx, Fy = c.CalculateFlux(k, ie, c.Q_Face)
					for n := 0; n < 4; n++ {
						normalFlux[i][n] = normal[0]*Fx[n] + normal[1]*Fy[n]
					}
				}
				c.SetNormalFluxOnRTEdge(k, edgeNumber, normalFlux, e.IInII[0])
			}
		case 2: // Handle edges with two connected tris - shared faces
			var (
				kL, kR                   = int(e.ConnectedTris[0]), int(e.ConnectedTris[1])
				edgeNumberL, edgeNumberR = int(e.ConnectedTriEdgeNumber[0]), int(e.ConnectedTriEdgeNumber[1])
				shiftL, shiftR           = edgeNumberL * Nedge, edgeNumberR * Nedge
				fluxLeft, fluxRight      [2][4]float64
			)
			switch c.FluxCalcAlgo {
			case FLUX_Average:
				averageFluxN := func(f1, f2 [2][4]float64, normal [2]float64) (fnorm [4]float64, fnormR [4]float64) {
					var (
						fave [2][4]float64
					)
					for n := 0; n < 4; n++ {
						for ii := 0; ii < 2; ii++ {
							fave[ii][n] = 0.5 * (f1[ii][n] + f2[ii][n])
						}
						fnorm[n] = normal[0]*fave[0][n] + normal[1]*fave[1][n]
						fnormR[n] = -fnorm[n]
					}
					return
				}
				normal, _ := c.getEdgeNormal(0, e, en)
				for i := 0; i < Nedge; i++ {
					iL := i + shiftL
					iR := Nedge - 1 - i + shiftR // Shared edges run in reverse order relative to each other
					fluxLeft[0], fluxLeft[1] = c.CalculateFlux(kL, iL, c.Q_Face)
					fluxRight[0], fluxRight[1] = c.CalculateFlux(kR, iR, c.Q_Face) // Reverse the right edge to match
					normalFlux[i], normalFluxReversed[Nedge-1-i] = averageFluxN(fluxLeft, fluxRight, normal)
				}
			case FLUX_LaxFriedrichs:
				var (
					rhoL, uL, vL, pL, CL float64
					rhoR, uR, vR, pR, CR float64
					maxV                 float64
				)
				maxVF := func(u, v, p, rho, C float64) (vmax float64) {
					vmax = math.Sqrt(u*u+v*v) + C
					return
				}
				normal, _ := c.getEdgeNormal(0, e, en)
				for i := 0; i < Nedge; i++ {
					iL := i + shiftL
					iR := Nedge - 1 - i + shiftR // Shared edges run in reverse order relative to each other
					rhoL, uL, vL, pL, CL, _ = c.GetState(kL, iL, c.Q_Face)
					rhoR, uR, vR, pR, CR, _ = c.GetState(kR, iR, c.Q_Face)
					fluxLeft[0], fluxLeft[1] = c.CalculateFlux(kL, iL, c.Q_Face)
					fluxRight[0], fluxRight[1] = c.CalculateFlux(kR, iR, c.Q_Face) // Reverse the right edge to match
					maxV = math.Max(maxVF(uL, vL, pL, rhoL, CL), maxVF(uR, vR, pR, rhoR, CR))
					indL, indR := kL+iL*Kmax, kR+iR*Kmax
					for n := 0; n < 4; n++ {
						normalFlux[i][n] = 0.5 * (normal[0]*(fluxLeft[0][n]+fluxRight[0][n]) + normal[1]*(fluxLeft[1][n]+fluxRight[1][n]))
						qD := c.Q_Face[n].Data()
						normalFlux[i][n] += 0.5 * maxV * (qD[indL] - qD[indR])
						normalFluxReversed[Nedge-1-i][n] = -normalFlux[i][n]
					}
				}
			case FLUX_Roe:
				var (
					rhoL, uL, vL, pL, EL float64
					rhoR, uR, vR, pR, ER float64
					hL, hR               float64
					Gamma                = c.Gamma
					GM1                  = Gamma - 1
				)
				normal, _ := c.getEdgeNormal(0, e, en)
				rotateMomentum := func(k, i int) {
					var (
						umD, vmD = c.Q_Face[1].Data(), c.Q_Face[2].Data()
					)
					ind := k + i*Kmax
					um, vm := umD[ind], vmD[ind]
					umD[ind] = um*normal[0] + vm*normal[1]
					vmD[ind] = -um*normal[1] + vm*normal[0]
				}
				for i := 0; i < Nedge; i++ {
					iL := i + shiftL
					iR := Nedge - 1 - i + shiftR // Shared edges run in reverse order relative to each other
					// Rotate the momentum into face normal coordinates before calculating fluxes
					rotateMomentum(kL, iL)
					rotateMomentum(kR, iR)
					rhoL, uL, vL, pL, _, EL = c.GetState(kL, iL, c.Q_Face)
					rhoR, uR, vR, pR, _, ER = c.GetState(kR, iR, c.Q_Face)
					fluxLeft[0], _ = c.CalculateFlux(kL, iL, c.Q_Face)
					fluxRight[0], _ = c.CalculateFlux(kR, iR, c.Q_Face) // Reverse the right edge to match
					/*
					   HM = (EnerM+pM).dd(rhoM);  HP = (EnerP+pP).dd(rhoP);
					*/
					// Enthalpy
					hL, hR = (EL+pL)/rhoL, (ER+pR)/rhoR
					/*
						rhoMs = sqrt(rhoM); rhoPs = sqrt(rhoP);
						rhoMsPs = rhoMs + rhoPs;

						rho = rhoMs.dm(rhoPs);
						u   = (rhoMs.dm(uM) + rhoPs.dm(uP)).dd(rhoMsPs);
						v   = (rhoMs.dm(vM) + rhoPs.dm(vP)).dd(rhoMsPs);
						H   = (rhoMs.dm(HM) + rhoPs.dm(HP)).dd(rhoMsPs);
						c2 = gm1 * (H - 0.5*(sqr(u)+sqr(v)));
						c = sqrt(c2);
					*/
					// Compute Roe average variables
					rhoLs, rhoRs := math.Sqrt(rhoL), math.Sqrt(rhoR)
					rhoLsRs := rhoLs + rhoRs

					rho := rhoLs * rhoRs
					u := (rhoLs*uL + rhoRs*uR) / rhoLsRs
					v := (rhoLs*vL + rhoRs*vR) / rhoLsRs
					h := (rhoLs*hL + rhoRs*hR) / rhoLsRs
					c2 := GM1 * (h - 0.5*(u*u+v*v))
					c := math.Sqrt(c2)
					/*
					   dW1 = -0.5*(rho.dm(uP-uM)).dd(c) + 0.5*(pP-pM).dd(c2);
					   dW2 = (rhoP-rhoM) - (pP-pM).dd(c2);
					   dW3 = rho.dm(vP-vM);
					   dW4 = 0.5*(rho.dm(uP-uM)).dd(c) + 0.5*(pP-pM).dd(c2);

					   dW1 = abs(u-c).dm(dW1);
					   dW2 = abs(u  ).dm(dW2);
					   dW3 = abs(u  ).dm(dW3);
					   dW4 = abs(u+c).dm(dW4);
					*/
					// Riemann fluxes
					dW1 := -0.5*(rho*(uR-uL))/c + 0.5*(pR-pL)/c2
					dW2 := (rhoR - rhoL) - (pR-pL)/c2
					dW3 := rho * (vR - vL)
					dW4 := 0.5*(rho*(uR-uL))/c + 0.5*(pR-pL)/c2
					dW1 = math.Abs(u-c) * dW1
					dW2 = math.Abs(u) * dW2
					dW3 = math.Abs(u) * dW3
					dW4 = math.Abs(u+c) * dW4
					/*
					   DMat fx = (fxQP+fxQM)/2.0;
					   fx(All,1) -= (dW1               + dW2                                   + dW4              )/2.0;
					   fx(All,2) -= (dW1.dm(u-c)       + dW2.dm(u)                             + dW4.dm(u+c)      )/2.0;
					   fx(All,3) -= (dW1.dm(v)         + dW2.dm(v)                 + dW3       + dW4.dm(v)        )/2.0;
					   fx(All,4) -= (dW1.dm(H-u.dm(c)) + dW2.dm(sqr(u)+sqr(v))/2.0 + dW3.dm(v) + dW4.dm(H+u.dm(c)))/2.0;
					*/
					// Form Roe Fluxes
					for n := 0; n < 4; n++ {
						normalFlux[i][n] = 0.5 * (fluxLeft[0][n] + fluxRight[0][n]) // Ave of normal component of flux
					}
					normalFlux[i][0] -= 0.5 * (dW1 + dW2 + dW4)
					normalFlux[i][1] -= 0.5 * (dW1*(u-c) + dW2*u + dW4*(u+c))
					normalFlux[i][2] -= 0.5 * (dW1*v + dW2*v + dW3 + dW4*v)
					normalFlux[i][3] -= 0.5 * (dW1*(h-u*c) + 0.5*dW2*(u*u+v*v) + dW3*v + dW4*(h+u*c))
					/*
					   flux = fx;    fx2.borrow(Ngf, fx.pCol(2)); fx3.borrow(Ngf, fx.pCol(3));
					   flux(All,2) = lnx.dm(fx2) - lny.dm(fx3);
					   flux(All,3) = lny.dm(fx2) + lnx.dm(fx3);
					*/
					// rotate back to Cartesian
					normalFlux[i][1], normalFlux[i][2] = normal[0]*normalFlux[i][1]-normal[1]*normalFlux[i][2], normal[1]*normalFlux[i][1]+normal[0]*normalFlux[i][2]
					// Project onto normal
					for n := 0; n < 4; n++ {
						normalFluxReversed[Nedge-1-i][n] = -normalFlux[i][n]
					}
				}
			}
			c.SetNormalFluxOnRTEdge(kL, edgeNumberL, normalFlux, e.IInII[0])
			c.SetNormalFluxOnRTEdge(kR, edgeNumberR, normalFluxReversed, e.IInII[1])
		}
	}
	return
}

func (c *Euler) getEdgeNormal(conn int, e *DG2D.Edge, en types.EdgeKey) (normal, scaledNormal [2]float64) {
	var (
		dfr = c.dfr
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
	revDir := bool(e.ConnectedTriDirection[conn])
	x1, x2 := DG2D.GetEdgeCoordinates(en, revDir, dfr.VX, dfr.VY)
	dx := [2]float64{x2[0] - x1[0], x2[1] - x1[1]}
	normal = normalize([2]float64{-dx[1], dx[0]})
	scaledNormal[0] = normal[0] * e.IInII[conn]
	scaledNormal[1] = normal[1] * e.IInII[conn]
	return
}

func (c *Euler) ProjectFluxToEdge(edgeFlux [][2][4]float64, e *DG2D.Edge, en types.EdgeKey) {
	/*
				Projects a 2D flux: [F, G] onto the face normal, then multiplies by the element/edge rqtio of normals, ||n||
		 		And places the scaled and projected normal flux into the degree of freedom F_RT_DOF
	*/
	var (
		dfr        = c.dfr
		Nedge      = dfr.FluxElement.Nedge
		Nint       = dfr.FluxElement.Nint
		Kmax       = dfr.K
		conn       = 0
		k          = int(e.ConnectedTris[conn])
		edgeNumber = int(e.ConnectedTriEdgeNumber[conn])
		shift      = edgeNumber * Nedge
	)
	// Get scaling factor ||n|| for each edge, multiplied by untransformed normals
	_, scaledNormal := c.getEdgeNormal(conn, e, en)
	for i := 0; i < Nedge; i++ {
		iL := i + shift
		// Project the flux onto the scaled scaledNormal
		for n := 0; n < 4; n++ {
			// Place normed/scaled flux into the RT element space
			rtD := c.F_RT_DOF[n].Data()
			ind := k + (2*Nint+iL)*Kmax
			rtD[ind] = scaledNormal[0]*edgeFlux[i][0][n] + scaledNormal[1]*edgeFlux[i][1][n]
		}
	}
}

func (c *Euler) SetNormalFluxOnRTEdge(k, edgeNumber int, edgeNormalFlux [][4]float64, IInII float64) {
	/*
		Takes the normal flux (aka "projected flux") multiplies by the ||n|| ratio of edge normals and sets that value for
		the F_RT_DOF degree of freedom locations for this [k, edgenumber] group
	*/
	var (
		dfr   = c.dfr
		Nedge = dfr.FluxElement.Nedge
		Nint  = dfr.FluxElement.Nint
		Kmax  = dfr.K
		shift = edgeNumber * Nedge
	)
	// Get scaling factor ||n|| for each edge, multiplied by untransformed normals
	for n := 0; n < 4; n++ {
		rtD := c.F_RT_DOF[n].Data()
		for i := 0; i < Nedge; i++ {
			// Place normed/scaled flux into the RT element space
			ind := k + (2*Nint+i+shift)*Kmax
			rtD[ind] = edgeNormalFlux[i][n] * IInII
		}
	}
}

func (c *Euler) InitializeFS() {
	var (
		K  = c.dfr.K
		Np = c.dfr.SolutionElement.Np
	)
	c.Q[0] = utils.NewMatrix(Np, K).AddScalar(c.Qinf[0])
	c.Q[1] = utils.NewMatrix(Np, K).AddScalar(c.Qinf[1])
	c.Q[2] = utils.NewMatrix(Np, K).AddScalar(c.Qinf[2])
	c.Q[3] = utils.NewMatrix(Np, K).AddScalar(c.Qinf[3])
}

func (c *Euler) InitializeIVortex(X, Y utils.Matrix) (iv *isentropic_vortex.IVortex, Q [4]utils.Matrix) {
	var (
		Beta     = 5.
		X0, Y0   = 5., 0.
		Gamma    = 1.4
		Np, Kmax = X.Dims()
	)
	for n := 0; n < 4; n++ {
		Q[n] = utils.NewMatrix(Np, Kmax)
	}
	iv = isentropic_vortex.NewIVortex(Beta, X0, Y0, Gamma)
	for ii := 0; ii < Np*Kmax; ii++ {
		x, y := X.Data()[ii], Y.Data()[ii]
		rho, rhoU, rhoV, E := iv.GetStateC(0, x, y)
		Q[0].Data()[ii] = rho
		Q[1].Data()[ii] = rhoU
		Q[2].Data()[ii] = rhoV
		Q[3].Data()[ii] = E
	}
	return
}

func (c *Euler) InitializeMemory() {
	var (
		K          = c.dfr.K
		Nedge      = c.dfr.FluxElement.Nedge
		NpFlux     = c.dfr.FluxElement.Np
		NpSolution = c.dfr.SolutionElement.Np
	)
	for n := 0; n < 4; n++ {
		c.Q_Face[n] = utils.NewMatrix(Nedge*3, K)
		c.F_RT_DOF[n] = utils.NewMatrix(NpFlux, K)
		c.DT = utils.NewMatrix(NpSolution, K) // Local time step
	}
}

func (c *Euler) CalculateFluxTransformed(k, i int, Q [4]utils.Matrix) (Fr, Fs [4]float64) {
	var (
		Jdet = c.dfr.Jdet.At(k, 0)
		Jinv = c.dfr.Jinv.Data()[4*k : 4*(k+1)]
	)
	Fx, Fy := c.CalculateFlux(k, i, Q)
	for n := 0; n < 4; n++ {
		Fr[n] = Jdet * (Jinv[0]*Fx[n] + Jinv[1]*Fy[n])
		Fs[n] = Jdet * (Jinv[2]*Fx[n] + Jinv[3]*Fy[n])
	}
	return
}

func (c *Euler) CalculateFlux(k, i int, Q [4]utils.Matrix) (Fx, Fy [4]float64) {
	// From https://www.theoretical-physics.net/dev/fluid-dynamics/euler.html
	var (
		q0D, q1D, q2D, q3D = Q[0].Data(), Q[1].Data(), Q[2].Data(), Q[3].Data()
		Kmax               = c.dfr.K
		ind                = k + i*Kmax
		rho, rhoU, rhoV, E = q0D[ind], q1D[ind], q2D[ind], q3D[ind]
	)
	Fx, Fy = c.FluxCalcMock(c.Gamma, rho, rhoU, rhoV, E)
	return
}

func FluxCalc(Gamma, rho, rhoU, rhoV, E float64) (Fx, Fy [4]float64) {
	var (
		GM1 = Gamma - 1.
	)
	u := rhoU / rho
	v := rhoV / rho
	u2 := u*u + v*v
	q := 0.5 * rho * u2
	p := GM1 * (E - q)
	Fx, Fy =
		[4]float64{rhoU, rhoU*u + p, rhoU * v, u * (E + p)},
		[4]float64{rhoV, rhoV * u, rhoV*v + p, v * (E + p)}
	return
}

func (c *Euler) RiemannBC(k, i int, QInf [4]float64, normal [2]float64) (Q [4]float64) {
	/*
			Use Riemann invariants along characteristic lines to calculate 1D flow properties normal to boundary
		Rinf = VnormInf - 2 * Cinf / (Gamma -1)
		Rint = VnormInt + 2 * Cint / (Gamma -1)
		Vn = 0.5 * (Rint + Rinf)
		C = 0.25 * (Gamma -1) *(Rint - Rinf)
		Then, project entropy and lateral velocity from the interior, calculate the other primitive variables from that
		P/(rho^Gamma) = constant
	*/
	var (
		rhoInt, uInt, vInt, pInt, CInt, _ = c.GetState(k, i, c.Q_Face)
		Gamma                             = c.Gamma
		GM1                               = Gamma - 1.
		OOGM1                             = 1. / GM1
		rhoInf, uInf, vInf, pInf, CInf, _ = GetPrimitiveVariables(Gamma, QInf)
		Vtang, Beta                       float64
		tangent                           = [2]float64{-normal[1], normal[0]}
	)
	VnormInt := normal[0]*uInt + normal[1]*vInt
	VnormInf := normal[0]*uInf + normal[1]*vInf
	Rinf := VnormInf - 2.*CInf*OOGM1
	Rint := VnormInt + 2.*CInt*OOGM1
	Vnorm := 0.5 * (Rint + Rinf)
	C := 0.25 * GM1 * (Rint - Rinf)
	//fmt.Printf("normal = %8.5f,%8.5f\n", normal[0], normal[1])
	switch {
	case VnormInt < 0: // Inflow, entropy and tangent velocity from Qinf
		Vtang = tangent[0]*uInf + tangent[1]*vInf
		Beta = pInf / math.Pow(rhoInf, Gamma)
	case VnormInt >= 0: // Outflow, entropy and tangent velocity from Qint
		Vtang = tangent[0]*uInt + tangent[1]*vInt
		Beta = pInt / math.Pow(rhoInt, Gamma)
	}
	u := Vnorm*normal[0] + Vtang*tangent[0]
	v := Vnorm*normal[1] + Vtang*tangent[1]
	//fmt.Printf("uInt,vInt=%8.5f,%8.5f u,v=%8.5f,%8.5f\n", uInt, vInt, u, v)
	rho := math.Pow(C*C/(Gamma*Beta), OOGM1)
	p := Beta * math.Pow(rho, Gamma)
	Q[0] = rho
	Q[1] = rho * u
	Q[2] = rho * v
	Q[3] = p*OOGM1 + 0.5*rho*(u*u+v*v)
	return
}

/************ Graphics **************
 */
type PlotField uint8

func (pm PlotField) String() string {
	strings := []string{"Density", "XMomentum", "YMomentum", "Energy", "Mach"}
	return strings[int(pm)]
}

const (
	Density PlotField = iota
	XMomentum
	YMomentum
	Energy
	Mach
)

type PlotMeta struct {
	Plot                   bool
	Scale                  float64
	TranslateX, TranslateY float64
	Field                  PlotField
	FieldMinP, FieldMaxP   *float64 // nil if no forced min, max
	FrameTime              time.Duration
	StepsBeforePlot        int
	LineType               chart2d.LineType
}

/************** Parallelism */

func (c *Euler) split1D(iMax, threadNum int) (istart, iend int) {
	var (
		Npart = iMax / c.ParallelDegree
	)
	istart = threadNum * Npart
	iend = istart + Npart
	if threadNum == c.ParallelDegree-1 {
		iend = iMax
	}
	return
}

package DG2D

import (
	"fmt"
	"math"
	"strconv"
	"testing"
	"time"

	"github.com/notargets/gocfd/utils"

	utils2 "github.com/notargets/avs/utils"

	"github.com/notargets/avs/chart2d"
	graphics2D "github.com/notargets/avs/geometry"

	"github.com/stretchr/testify/assert"
)

func TestLagrangePolynomial(t *testing.T) {
	numSamples := 500
	rd := make([]float64, numSamples)
	xmin, xmax := -1., 1.
	fmin, fmax := -0.5, 1.25
	inc := (xmax - xmin) / float64(numSamples-1.)
	for i := 0; i < numSamples; i++ {
		rd[i] = xmin + float64(i)*inc
	}
	SamplesR := utils.NewVector(numSamples, rd)
	f := make([]float64, SamplesR.Len())
	var plot bool
	plot = false
	if plot {
		chart := utils.NewLineChart(1920, 1080, xmin, xmax, fmin, fmax)
		// TODO: Make a pluggable basis underneath the RT (and Lagrange) elements - Lagrange, Hesthaven, Spectral?
		var delay time.Duration
		Nmax := 4
		lineInc := 2. / float64(Nmax-2)
		lineColor := -1. // colormap goes from -1,1
		for n := 1; n < Nmax; n++ {
			Np := n + 1
			R := utils.NewVector(Np)
			inc = (xmax - xmin) / float64(Np-1)
			for i := 0; i < Np; i++ {
				R.DataP[i] = xmin + float64(i)*inc
			}
			lp := NewLagrangeBasis1D(R.DataP)
			_ = lp
			for j := 0; j < R.Len(); j++ {
				for i, r := range SamplesR.DataP {
					//f[i] = DG1D.JacobiP(utils.NewVector(1, []float64{r}), 0, 0, j)[0]
					f[i] = lp.EvaluateBasisPolynomial([]float64{r}, j)[0]
				}
				if n == Nmax-1 && j == R.Len()-1 {
					delay = 120 * time.Second
				}
				name := "JacobiP[" + strconv.Itoa(n) + "," + strconv.Itoa(j) + "]"
				fmt.Printf("Chart Name: [%s], lineColor = %5.3f\n", name, lineColor)
				chart.Plot(delay, SamplesR.DataP, f, lineColor, name)
			}
			lineColor += lineInc
		}
	}
	{
		R := []float64{-1, -.5, 0, .5, 1.} // P = 4
		lp := NewLagrangeBasis1D(R)
		RI := []float64{-.9, -.75, -.2, .2, .75, .9} // 6 points
		im := lp.GetInterpolationMatrix(RI)
		Ni, Np := im.Dims()
		assert.Equal(t, 6, Ni)
		assert.Equal(t, 5, Np)
		F := utils.NewMatrix(Np, 1, []float64{10, 5, 1, 5, 10})
		FI := im.Mul(F)
		assert.InDeltaSlice(t, []float64{9.72640, 8.35938, 1.71840, 1.71840, 8.35938, 9.72640},
			FI.DataP, 0.00001, "err msg %s")
	}
}

type LagrangeBasis1D struct {
	P       int       // Order
	Np      int       // Dimension of basis = N+1
	Weights []float64 // Barycentric weights, one per basis polynomial
	Nodes   []float64 // Nodes at which basis is defined
}

func NewLagrangeBasis1D(r []float64) (lb *LagrangeBasis1D) {
	/*
		At a given order P, there are (P+1) basis polynomials representing that order
		To recover a basis polynomial we need to specifiy:
		`	P = The order of the basis
			j = The basis polynomial number within the basis
			R = The points used to define the basis, (P+1) dimension
	*/
	lb = &LagrangeBasis1D{
		P:       len(r) - 1,
		Np:      len(r),
		Weights: make([]float64, len(r)),
		Nodes:   r,
	}
	// Calculate the weight for each basis function j
	for j := 0; j < lb.Np; j++ {
		lb.Weights[j] = 1.
	}
	for j := 0; j < lb.Np; j++ {
		for i := 0; i < lb.Np; i++ {
			if i != j {
				lb.Weights[j] /= r[j] - r[i]
			}
		}
	}
	return
}

func (lb *LagrangeBasis1D) GetInterpolationMatrix(R []float64) (im utils.Matrix) {
	/*
			Provided function values at each of the P+1 nodes, interpolate a new function value at location r
			Note that the points in R are not necessarily the defining points of the basis, and are not necessarily at the
		    same points within F, the provided set of function values at the nodes of the basis
	*/
	var (
		fj = make([]float64, len(R)) // temporary storage for each basis function evaluation
	)
	im = utils.NewMatrix(len(R), lb.Np) // Rows are for evaluation points, columns for basis
	for j := 0; j < lb.Np; j++ {        // For each basis function
		fj = lb.EvaluateBasisPolynomial(R, j)
		for i, val := range fj {
			im.Set(i, j, val)
		}
	}
	return
}

func (lb *LagrangeBasis1D) Interpolate(R []float64, F []float64) (f []float64) {
	/*
			Provided function values at each of the P+1 nodes, interpolate a new function value at location r
			Note that the points in R are not necessarily the defining points of the basis, and are not necessarily at the
		    same points within F, the provided set of function values at the nodes of the basis
	*/
	var (
		fj = make([]float64, len(R)) // temporary storage for each basis function evaluation
	)
	for j := 0; j < lb.Np; j++ { // For each basis function
		fj = lb.EvaluateBasisPolynomial(R, j)
		for i := range R {
			f[i] += fj[i] * F[j]
		}
	}
	return
}

func (lb *LagrangeBasis1D) EvaluateBasisPolynomial(R []float64, j int) (f []float64) {
	/*
		This evaluates a single basis polynomial (the jth) within the basis for order P at all points in R
		Note that the points in R are not necessarily the defining points of the basis
	*/
	f = make([]float64, len(R))
	for i, r := range R {
		f[i] = lb.evaluateL(r) * lb.Weights[j]
		if math.Abs(r-lb.Nodes[j]) < 0.0000000001 {
			f[i] = 1.
		} else {
			f[i] /= (r - lb.Nodes[j])
		}
	}
	return
}

func (lb *LagrangeBasis1D) evaluateL(r float64) (f float64) {
	/*
		This is the polynomial term in the Barycentric version of the Lagrange polynomial basis
		It is not specific to the jth polynomial, but applies to all the individual basis polynomials
	*/
	f = 1.
	for _, rr := range lb.Nodes {
		f *= (r - rr)
	}
	return
}

func TestRTElement(t *testing.T) {
	{
		// Check term-wise orthogonal 2D polynomial basis
		N := 2
		R, S := NodesEpsilon(N - 1)
		ii, jj := 1, 1
		p := Simplex2DP(R, S, ii, jj)
		ddr, dds := GradSimplex2DP(R, S, ii, jj)
		Np := R.Len()
		pCheck, ddrCheck, ddsCheck := make([]float64, Np), make([]float64, Np), make([]float64, Np)
		for i, rVal := range R.DataP {
			sVal := S.DataP[i]
			ddrCheck[i], ddsCheck[i] = GradSimplex2DPTerm(rVal, sVal, ii, jj)
			pCheck[i] = Simplex2DPTerm(rVal, sVal, ii, jj)
		}
		assert.True(t, nearVec(pCheck, p, 0.000001))
		assert.True(t, nearVec(ddrCheck, ddr, 0.000001))
		assert.True(t, nearVec(ddsCheck, dds, 0.000001))
	}
	errorCheck := func(N int, div, divCheck []float64) (minInt, maxInt, minEdge, maxEdge float64) {
		var (
			Npm    = len(div)
			errors = make([]float64, Npm)
		)
		for i := 0; i < Npm; i++ {
			//var ddr, dds float64
			errors[i] = div[i] - divCheck[i]
		}
		minInt, maxInt = errors[0], errors[0]
		Nint := N * (N + 1) / 2
		minEdge, maxEdge = errors[Nint], errors[Nint]
		for i := 0; i < Nint; i++ {
			errAbs := math.Abs(errors[i])
			if minInt > errAbs {
				minInt = errAbs
			}
			if maxInt < errAbs {
				maxInt = errAbs
			}
		}
		for i := Nint; i < Npm; i++ {
			errAbs := math.Abs(errors[i])
			if minEdge > errAbs {
				minEdge = errAbs
			}
			if maxEdge < errAbs {
				maxEdge = errAbs
			}
		}
		fmt.Printf("Order = %d, ", N)
		fmt.Printf("Min, Max Int Err = %8.5f, %8.5f, Min, Max Edge Err = %8.5f, %8.5f\n", minInt, maxInt, minEdge, maxEdge)
		return
	}
	checkSolution := func(rt *RTElement, Order int) (s1, s2, divCheck []float64) {
		var (
			Np = rt.Np
		)
		s1, s2 = make([]float64, Np), make([]float64, Np)
		divCheck = make([]float64, Np)
		var ss1, ss2 float64
		for i := 0; i < Np; i++ {
			r := rt.R.DataP[i]
			s := rt.S.DataP[i]
			ccf := float64(Order)
			s1[i] = utils.POW(r, Order)
			s2[i] = utils.POW(s, Order)
			ss1, ss2 = ccf*utils.POW(r, Order-1), ccf*utils.POW(s, Order-1)
			divCheck[i] = ss1 + ss2
		}
		return
	}
	{
		N := 1
		R, S := NodesEpsilon(N - 1)
		rt := NewRTElement(N, R, S)
		fmt.Println(rt.V[0].Print("V0"))
		fmt.Println(rt.V[1].Print("V1"))
		fmt.Println(rt.Div.Print("Div"))
	}
	if true { // Check Divergence for polynomial vector fields of order < N against analytical solution
		Nend := 8
		for N := 1; N < Nend; N++ {
			R, S := NodesEpsilon(N - 1)
			rt := NewRTElement(N, R, S)
			for cOrder := 0; cOrder <= N; cOrder++ {
				fmt.Printf("Check Order = %d, ", cOrder)
				s1, s2, divCheck := checkSolution(rt, cOrder)
				sp := rt.ProjectFunctionOntoBasis(s1, s2)
				sm := utils.NewMatrix(rt.Np, 1, sp)
				divM := rt.Div.Mul(sm)
				//fmt.Println(divM.Print("divM"))
				minerrInt, maxerrInt, minerrEdge, maxerrEdge := errorCheck(N, divM.DataP, divCheck)
				assert.True(t, near(minerrInt, 0.0, 0.00001))
				assert.True(t, near(maxerrInt, 0.0, 0.00001))
				assert.True(t, near(minerrEdge, 0.0, 0.00001))
				assert.True(t, near(maxerrEdge, 0.0, 0.00001))
			}
		}
	}
	plot := false
	if plot {
		N := 7
		NRT := N + 1
		R, S := NodesEpsilon(N)
		//Nint := R.Len()
		rt := NewRTElement(NRT, R, S)
		s1, s2 := make([]float64, rt.R.Len()), make([]float64, rt.R.Len())
		for i := range rt.R.DataP {
			s1[i] = 1
			s2[i] = 1
		}
		if plot {
			chart := PlotTestTri(true)
			points := utils.ArraysToPoints(rt.R.DataP, rt.S.DataP)
			f := utils.ArraysTo2Vector(s1, s2, 0.1)
			_ = chart.AddVectors("test function", points, f, chart2d.Solid, utils.GetColor(utils.Green))
			utils.SleepFor(500000)
		}
	}
}

func PlotTestTri(plotGeom bool) (chart *chart2d.Chart2D) {
	var (
		points  []graphics2D.Point
		trimesh graphics2D.TriMesh
		K       = 1
	)

	points = make([]graphics2D.Point, 3)
	points[0].X[0], points[0].X[1] = -1, -1
	points[1].X[0], points[1].X[1] = 1, -1
	points[2].X[0], points[2].X[1] = -1, 1

	trimesh.Triangles = make([]graphics2D.Triangle, K)
	colorMap := utils2.NewColorMap(0, 1, 1)
	trimesh.Triangles[0].Nodes[0] = 0
	trimesh.Triangles[0].Nodes[1] = 1
	trimesh.Triangles[0].Nodes[2] = 2
	trimesh.Geometry = points
	box := graphics2D.NewBoundingBox(trimesh.GetGeometry())
	chart = chart2d.NewChart2D(1024, 1024, box.XMin[0], box.XMax[0], box.XMin[1], box.XMax[1])
	chart.AddColorMap(colorMap)
	go chart.Plot()

	if plotGeom {
		if err := chart.AddTriMesh("TriMesh", trimesh,
			chart2d.CrossGlyph, chart2d.Solid, utils.GetColor(utils.White)); err != nil {
			panic("unable to add graph series")
		}
	}
	return
}

/*
	checkSolutionM := func(rt *RTElement, Order int) (s1, s2, divCheck []float64) {
		var (
			Npm = rt.Npm
		)
		s1, s2 = make([]float64, Npm), make([]float64, Npm)
		divCheck = make([]float64, Npm)
		var ss1, ss2 float64
		for i := 0; i < Npm; i++ {
			r := rt.Rm.DataP[i]
			s := rt.Sm.DataP[i]
			ccf := float64(Order)
			s1[i] = utils.POW(r, Order)
			s2[i] = utils.POW(s, Order)
			ss1, ss2 = ccf*utils.POW(r, Order-1), ccf*utils.POW(s, Order-1)
			divCheck[i] = ss1 + ss2
		}
		return
	}
	if false { // Check Divergence for polynomial vector fields of order < N against analytical solution
		Nend := 2
		for N := 1; N < Nend; N++ {
			R, S := NodesEpsilon(N - 1)
			rt := NewRTElement(N, R, S)
			fmt.Printf("calculating Dr and Ds...\n")
			//Dr := rt.Vrm[0].Mul(rt.VmInv[0])
			//Ds := rt.Vsm[1].Mul(rt.VmInv[1])
			Dr := rt.Vrm[0].Mul(rt.VmInv[0])
			Ds := rt.Vsm[1].Mul(rt.VmInv[1])
			//Dr := rt.Vrm[0]
			//Ds := rt.Vsm[1]
			for cOrder := 0; cOrder <= N; cOrder++ {
				fmt.Printf("Check Order = %d, ", cOrder)
				s1, s2, divCheck := checkSolutionM(rt, cOrder)
				s1p, s2p := rt.ProjectFunctionOntoBasis2(s1, s2)
				s1m, s2m := utils.NewMatrix(rt.Npm, 1, s1p), utils.NewMatrix(rt.Npm, 1, s2p)
				div := Dr.Mul(s1m).Add(Ds.Mul(s2m)).DataP
				minerrInt, maxerrInt, minerrEdge, maxerrEdge := errorCheck(N, div, divCheck)
				assert.True(t, near(minerrInt, 0.0, 0.00001))
				assert.True(t, near(maxerrInt, 0.0, 0.00001))
				assert.True(t, near(minerrEdge, 0.0, 0.00001))
				assert.True(t, near(maxerrEdge, 0.0, 0.00001))
			}
		}
	}
*/

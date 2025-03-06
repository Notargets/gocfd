package DG2D

import (
	"fmt"
	"math"
	"sort"
	"testing"
	"time"

	"github.com/notargets/avs/geometry"

	"github.com/notargets/gocfd/geometry2D"

	"github.com/notargets/gocfd/model_problems/Euler2D/isentropic_vortex"

	"github.com/notargets/avs/chart2d"
	utils2 "github.com/notargets/avs/utils"

	"github.com/notargets/gocfd/utils"
)

type TestField uint8

const (
	NORMALSHOCKTESTM2 = TestField(iota)
	NORMALSHOCKTESTM5 = TestField(iota)
	NORMALSHOCKTESTM12
	FIXEDVORTEXTEST
	MOVINGVORTEXTEST
	RADIAL1TEST
	RADIAL2TEST
	RADIAL3TEST
	RADIAL4TEST
)

func (tf TestField) String() string {
	switch tf {
	case NORMALSHOCKTESTM2:
		return "NORMALSHOCKTESTM2"
	case NORMALSHOCKTESTM5:
		return "NORMALSHOCKTESTM5"
	case NORMALSHOCKTESTM12:
		return "NORMALSHOCKTESTM12"
	case FIXEDVORTEXTEST:
		return "FIXEDVORTEXTEST"
	case MOVINGVORTEXTEST:
		return "MOVINGVORTEXTEST"
	case RADIAL1TEST:
		return "RADIAL1TEST"
	case RADIAL2TEST:
		return "RADIAL2TEST"
	case RADIAL3TEST:
		return "RADIAL3TEST"
	case RADIAL4TEST:
		return "RADIAL4TEST"
	}
	return ""
}

func setTestField(X, Y []float64, tf TestField) (field []float64) {
	switch tf {
	case NORMALSHOCKTESTM12:
		field = setShockConditions(X, 1.2, 0)
	case NORMALSHOCKTESTM2:
		field = setShockConditions(X, 2, 0)
	case NORMALSHOCKTESTM5:
		field = setShockConditions(X, 5, 0)
	case FIXEDVORTEXTEST:
		iv := isentropic_vortex.NewIVortex(5, 0, 0, 1.4, 0)
		field = setIsoVortexConditions(X, Y, iv)
	case MOVINGVORTEXTEST:
		iv := isentropic_vortex.NewIVortex(5, 0, 0, 1.4)
		field = setIsoVortexConditions(X, Y, iv)
	case RADIAL1TEST:
		field = setRadial(X, Y, 1)
	case RADIAL2TEST:
		field = setRadial(X, Y, 2)
	case RADIAL3TEST:
		field = setRadial(X, Y, 3)
	case RADIAL4TEST:
		field = setRadial(X, Y, 4)
	}
	return
}

func convXYtoXandY(XY []float32) (X, Y []float64) {
	for i := 0; i < len(XY)/2; i++ {
		X = append(X, float64(XY[2*i]))
		Y = append(Y, float64(XY[2*i+1]))
	}
	return
}

func setRadial(X, Y []float64, order float64) (field []float64) {
	var (
		Np = len(X)
	)
	field = make([]float64, Np*4)
	var x, y float64
	for i := 0; i < Np; i++ {
		x, y = X[i], Y[i]
		for n := 0; n < 4; n++ {
			field[i+n*Np] = math.Pow(x, order) + math.Pow(y, order)
		}
	}
	return
}

func setShockConditions(X []float64, Mach, Alpha float64) (field []float64) {
	var (
		U1, U2 = ShockConditions(Mach, Alpha)
		Np     = len(X)
	)
	// fmt.Printf("Mach: %.2f\n", Mach)
	// fmt.Printf("U1: ")
	// for n := 0; n < 4; n++ {
	// 	fmt.Printf(" %.2f ", U1[n])
	// }
	// fmt.Printf("\n")
	// fmt.Printf("U2: ")
	// for n := 0; n < 4; n++ {
	// 	fmt.Printf(" %.2f ", U2[n])
	// }
	// fmt.Printf("\n")
	field = make([]float64, Np*4)
	var x float64
	for i := 0; i < Np; i++ {
		x = X[i]
		for n := 0; n < 4; n++ {
			if x < 0 {
				field[i+n*Np] = U1[n]
			} else {
				field[i+n*Np] = U2[n]
			}
		}
	}
	return
}

func setIsoVortexConditions(X, Y []float64,
	iv *isentropic_vortex.IVortex) (field []float64) {
	var (
		Np = len(X)
	)
	field = make([]float64, Np*4)
	var x, y float64
	for i := 0; i < Np; i++ {
		x, y = X[i], Y[i]
		rho, rhoU, rhoV, rhoE := iv.GetStateC(0, x, y)
		field[i+0*Np] = rho
		field[i+1*Np] = rhoU
		field[i+2*Np] = rhoV
		field[i+3*Np] = rhoE
	}
	return
}

func TestPlotVariousFields(t *testing.T) {
	var (
		N = 2
	)
	if !testing.Verbose() {
		return
	}
	angle := 210.
	// angle := 0.
	dfr := CreateEquiTriMesh(N, angle)
	gm := CreateGraphMesh(dfr)
	X, Y := convXYtoXandY(gm.XY)
	field := setTestField(X, Y, FIXEDVORTEXTEST)
	// field := setTestField(X, Y, NORMALSHOCKTESTM2)
	// field := setTestField(X, Y, RADIAL2TEST)
	n := 0 // Density
	n = 1  // U momentum
	n = 2  // V momentum
	n = 3  // rho*E energy
	n = 0  // Density
	plotField(field, n, gm)
}

func TestInterpolationVariousFields(t *testing.T) {
	var (
		NMin  = 1
		NMax  = 7
		Nu, p = 0.2, 3.
	)
	for N := NMin; N <= NMax; N++ {
		fmt.Printf("ORDER: %d Element Test\n-----------------------\n", N)
		// angle := 210.
		angle := 0.
		dfr := CreateEquiTriMesh(N, angle)
		rt := dfr.FluxElement
		RFlux := utils.NewVector(rt.NpEdge*3, rt.GetEdgeLocations(rt.R.DataP)) // For the Interpolation matrix across three edges
		SFlux := utils.NewVector(rt.NpEdge*3, rt.GetEdgeLocations(rt.S.DataP)) // For the Interpolation matrix across three edges
		EdgeInterpMod := dfr.SolutionElement.JB2D.GetModInterpMatrix(RFlux,
			SFlux, Nu, p)
		for _, tf := range []TestField{NORMALSHOCKTESTM12, NORMALSHOCKTESTM2,
			NORMALSHOCKTESTM5, FIXEDVORTEXTEST, RADIAL1TEST, RADIAL2TEST,
			RADIAL3TEST, RADIAL4TEST} {
			// for _, tf := range []TestField{NORMALSHOCKTESTM12} {
			Np := dfr.SolutionElement.Np
			X, Y := dfr.SolutionX.DataP, dfr.SolutionY.DataP
			QSol := QFromField(setTestField(X, Y, tf), Np)
			fmt.Printf("%s Interpolation\n", tf.String())
			RMSErr, Mean, trueEdges, interpEdges :=
				GetInterpolationAccuracy(dfr, QSol, dfr.FluxEdgeInterp, tf)
			printRMSError(RMSErr, Mean, N)

			fmt.Printf("Modulated Field\n")
			QSolMod := ModulateInternalField(dfr, QSol, Nu, p, 5)
			RMSErr, Mean, trueEdges, interpEdges =
				GetInterpolationAccuracy(dfr, QSolMod, dfr.FluxEdgeInterp, tf)
			printRMSError(RMSErr, Mean, N)

			fmt.Printf("Modulated Field and Modulated Interpolation\n")
			RMSErr, Mean, trueEdges, interpEdges =
				GetInterpolationAccuracy(dfr, QSolMod, EdgeInterpMod, tf)
			printRMSError(RMSErr, Mean, N)
			_, _ = trueEdges, interpEdges
			fmt.Printf("---------------------\n\n\n")
			// Coeffs := GetRTCoefficients(dfr, tf)
			// QSol.Print("QSol Orig")
			// for _, iterCount := range []int{10, 100, 1000} {
			// 	QSolMod := ModulateInternalField(dfr, QSol, Nu, p, iterCount)
			// 	QSolMod.Print("QSolMod." + strconv.Itoa(iterCount))
			// }
		}
	}
}

func (jb2d *JacobiBasis2D) GetModInterpMatrix(R, S utils.Vector,
	Nu, p float64) (Interp utils.Matrix) {
	/*
		Uses Jacobi polynomials as the basis function

		Compose a matrix of interpolating polynomials where each row represents one [r,s] location to be interpolated
		This matrix can then be multiplied by a single vector of function values at the polynomial nodes to produce a
		vector of interpolated values, one for each interpolation location
	*/
	var (
		N  = jb2d.P
		Np = jb2d.Np
	)
	CoefMods := ModalFilter2D(Nu, p, N)
	// First compute polynomial terms, used by all polynomials
	polyTerms := make([]float64, R.Len()*Np)
	var sk int
	for ii, r := range R.DataP {
		s := S.DataP[ii]
		var sk2 int
		for i := 0; i <= N; i++ {
			for j := 0; j <= (N - i); j++ {
				polyTerms[sk] = CoefMods[sk2] * jb2d.PolynomialTerm(r, s, i, j)
				sk++
				sk2++
			}
		}
	}
	ptV := utils.NewMatrix(R.Len(), Np, polyTerms).Transpose()
	Interp = jb2d.Vinv.Transpose().Mul(ptV).Transpose()
	return
}

func ModulateInternalField(dfr *DFR2D, QSol utils.Matrix, Nu,
	p float64, iterCount int) (QSolMod utils.Matrix) {
	var (
		Np = dfr.SolutionElement.Np
	)
	QSolMod = utils.NewMatrix(Np, 4)
	CoeffModifier := ModalFilter2D(Nu, p, dfr.N)
	for iter := 0; iter < iterCount; iter++ {
		for n := 0; n < 4; n++ {
			Coeffs := dfr.SolutionElement.JB2D.Vinv.Mul(QSol.Col(n).ToMatrix())
			for i := 0; i < Np; i++ {
				Coeffs.DataP[i] *= CoeffModifier[i]
			}
			QSolModE := dfr.SolutionElement.JB2D.V.Mul(Coeffs)
			for i := 0; i < Np; i++ {
				QSolMod.Set(i, n, QSolModE.DataP[i])
			}
		}
		QSol = QSolMod.Copy(QSol)
	}
	return
}

func ModalFilter2D(Nu, p float64, P int) (CoeffModifier []float64) {
	// Tunables:
	// - Nu is a user-defined dissipation strength (tunable),
	// 0.1 to 0.5	Higher = stronger overall dissipation (aggressiveness)
	// - p controls sharpness (typically 1 to 4 — sharper if you only want to hit the highest modes).
	// 1 to 4	Higher = sharper cutoff (more selective to highest modes)
	var (
		degree = ModalDegree2D(P)
		Np     = len(degree)
	)
	CoeffModifier = make([]float64, Np)
	for i := 0; i < Np; i++ {
		CoeffModifier[i] = 1. - Nu*(math.Pow(degree[i]/float64(P), 2.*p))
	}
	return
}

func ModalDegree2D(P int) (degree []float64) {
	var (
		Np = (P + 1) * (P + 2) / 2
	)
	degree = make([]float64, Np)
	var sk int
	for i := 0; i <= P; i++ {
		for j := 0; j <= P-i; j++ {
			degree[sk] = math.Max(float64(i), float64(j))
			sk++
		}
	}
	return
}

func GetRTCoefficients(dfr *DFR2D, tf TestField) (Coeffs [4]utils.Matrix) {
	var (
		Np            = dfr.FluxElement.Np
		X, Y          = dfr.FluxX.DataP, dfr.FluxY.DataP
		QFlux         = QFromField(setTestField(X, Y, tf), Np)
		F, G          = GetFluxVectors(QFlux)
		ProjectedFlux = utils.NewMatrix(Np, 4)
	)
	for n := 0; n < 4; n++ {
		b, e := n*Np, (n+1)*Np
		dfr.FluxElement.ProjectFunctionOntoDOF(F.DataP[b:e], G.DataP[b:e],
			ProjectedFlux.DataP[b:e])
		Coeffs[n] = dfr.FluxElement.VInv.Mul(ProjectedFlux.Col(n).ToMatrix())
	}
	return
}

func GetFluxVectors(QFlux utils.Matrix) (F, G utils.Matrix) {
	var (
		gamma = 1.4
		Np, _ = QFlux.Dims()
	)
	F, G = utils.NewMatrix(Np, 4), utils.NewMatrix(Np, 4)
	P := func(rho, u, v, E float64) (p float64) {
		// p=(γ−1)(E−(ρ/2)(u^2+v^2))
		p = (gamma - 1.) * (E - (rho/2.)*(u*u+v*v))
		return
	}
	for i := 0; i < Np; i++ {
		rho, rhou, rhov, E := QFlux.At(i, 0), QFlux.At(i, 1), QFlux.At(i, 2),
			QFlux.At(i, 3)
		u, v := rhou/rho, rhov/rho
		p := P(rho, u, v, E)
		F.Set(i, 0, rhou)
		G.Set(i, 0, rhov)

		F.Set(i, 1, u*rhou+p)
		G.Set(i, 1, u*rhov)

		F.Set(i, 2, v*rhou)
		G.Set(i, 2, v*rhov+p)

		F.Set(i, 3, u*(E+p))
		G.Set(i, 3, v*(E+p))
	}
	return
}

func GetInterpolationAccuracy(dfr *DFR2D, QSol, EdgeInterpolation utils.Matrix,
	tf TestField) (RMSErr, Mean [4]float64, trueEdges, interpEdges [4][]float64) {
	var (
		NpInt          = dfr.FluxElement.NpInt
		edgeX          = dfr.FluxX.DataP[2*NpInt:]
		edgeY          = dfr.FluxY.DataP[2*NpInt:]
		trueEdgesField = setTestField(edgeX, edgeY, tf)
		fieldLen       = len(trueEdgesField) / 4
	)
	interpEdges, RMSErr = evaluateInterpolation(QSol, EdgeInterpolation,
		trueEdgesField)
	for n := 0; n < 4; n++ {
		trueEdges[n] = make([]float64, fieldLen)
		for i := 0; i < fieldLen; i++ {
			trueEdges[n][i] = trueEdgesField[i+n*fieldLen]
		}
		for _, f := range trueEdges[n] {
			Mean[n] += f
		}
		Mean[n] /= float64(len(trueEdges[n]))
		if math.Abs(Mean[n]) < 0.0001 {
			Mean[n] = 1.
		}
	}
	return
}

func printRMSError(RMSErr, Mean [4]float64, N int) {
	fmt.Printf("Order:%d RMS: ", N)
	for nn := 0; nn < 4; nn++ {
		fmt.Printf(" %.2f,", RMSErr[nn])
	}
	fmt.Printf("\tRMS%%: ")
	for nn := 0; nn < 4; nn++ {
		fmt.Printf(" %.2f%%,", 100.*RMSErr[nn]/Mean[nn])
	}
	fmt.Printf("\n")
}

func evaluateInterpolation(QSol utils.Matrix, EdgeInterpolation utils.Matrix,
	trueVals []float64,
	printO ...interface{}) (edges [4][]float64, RMSErr [4]float64) {
	var (
		MinMaxSol  [4][2]float64
		MinMaxInt  [4][2]float64
		MinMaxTrue [4][2]float64
		MEAN       [4]float64
		lenField   = len(trueVals) / 4
	)
	minmax := func(f []float64) (min, max float64) {
		min, max = f[0], f[0]
		for _, ff := range f[1:] {
			if min > ff {
				min = ff
			}
			if max < ff {
				max = ff
			}
		}
		return
	}
	for n := 0; n < 4; n++ {
		tField := trueVals[n*lenField : (n+1)*lenField]
		// edges[n] = dfr.FluxEdgeInterp.Mul(QSol.Col(n).ToMatrix()).DataP
		edges[n] = EdgeInterpolation.Mul(QSol.Col(n).ToMatrix()).DataP
		MinMaxSol[n][0] = QSol.Col(n).Min()
		MinMaxSol[n][1] = QSol.Col(n).Max()
		MinMaxInt[n][0], MinMaxInt[n][1] = minmax(edges[n])
		MinMaxTrue[n][0], MinMaxTrue[n][1] = minmax(tField)
		for i := 0; i < lenField; i++ {
			e := edges[n][i] - tField[i]
			RMSErr[n] += e * e
			MEAN[n] += tField[i]
		}
		RMSErr[n] /= float64(lenField)
		RMSErr[n] = math.Sqrt(RMSErr[n])
		MEAN[n] /= float64(lenField)
		if math.Abs(MEAN[n]) < 0.001 {
			MEAN[n] = 1
		}
	}
	printField := func(field []float64, n int, label string,
		minmax [2]float64, rms ...float64) {
		label2 := label[0:1]
		fmt.Printf("%s[%d]\tQMin:%+.1f\tQMax:%+.1f\t"+
			"%sMin:%+.1f\t%sMax:%+.1f\t",
			label, n, MinMaxSol[n][0], MinMaxSol[n][1],
			label2, minmax[0], label2, minmax[1])
		for _, f := range field {
			if math.Abs(f) < 0.0001 {
				f = 0.
			}
			fmt.Printf(" %.1f ", f)
		}
		if len(rms) > 0 {
			fmt.Printf(" RMSErr = %.2f, %%RMS/MEAN = %.2f%%\n",
				RMSErr[n], 100*RMSErr[n]/MEAN[n])
		} else {
			fmt.Printf("\n")
		}
	}
	if len(printO) != 0 {
		for n := 0; n < 4; n++ {
			printField(edges[n], n, "Interp", MinMaxInt[n])
			printField(trueVals[n*lenField:(n+1)*lenField], n, "True  ",
				MinMaxTrue[n], RMSErr[n])
		}
	}
	return
}

func QFromField(field []float64, Np int) (Q utils.Matrix) {
	Q = utils.NewMatrix(Np, 4)
	for n := 0; n < 4; n++ {
		for i := 0; i < Np; i++ {
			Q.Set(i, n, field[i+Np*n])
		}
	}
	return
}

func plotField(field []float64, fieldNum int, gm geometry.TriMesh,
	xMM ...float64) {
	var xMin, xMax, yMin, yMax float32

	if len(xMM) == 4 {
		xMin, xMax = float32(xMM[0]), float32(xMM[1])
		yMin, yMax = float32(xMM[2]), float32(xMM[3])
	} else {
		xMin, xMax = -1, 1
		yMin, yMax = -1, 1
	}
	ch := chart2d.NewChart2D(xMin, xMax, yMin, yMax,
		1024, 1024, utils2.WHITE, utils2.BLACK)
	// Create a vector field including the three vertices
	var pField []float32
	vertsLen := len(field) / 4
	var fMin, fMax float32
	fMin, fMax = math.MaxFloat32, -math.MaxFloat32
	n := 0 // Density
	n = 2  // V momentum
	n = 3  // rho*E energy
	n = 0  // Density
	n = 1  // U momentum
	n = 3  // rho*E energy
	n = 0  // Density
	for i := 0; i < vertsLen; i++ {
		pField = append(pField, float32(field[i+n*vertsLen]))
		// fmt.Printf("F[%d] = %.2f\n", i, pField[i])
		if pField[i] < fMin {
			fMin = pField[i]
		}
		if pField[i] > fMax {
			fMax = pField[i]
		}
	}
	vs := geometry.VertexScalar{
		TMesh:       &gm,
		FieldValues: pField,
	}
	_ = vs
	_ = ch
	// ch.AddContourVertexScalar(&vs, fMin, fMax, 100)
	ch.AddShadedVertexScalar(&vs, fMin, fMax)
	ch.AddTriMesh(gm)
	for {
	}
}

func CreateGraphMesh(dfr *DFR2D) (gm geometry.TriMesh) {
	var (
		NpInt        = dfr.FluxElement.NpInt
		NpEdge       = dfr.FluxElement.NpEdge
		TriNp        = 3 + 3*NpEdge
		EdgeX, EdgeY = make([]float64, TriNp), make([]float64, TriNp)
	)
	// Compose the edges of the triangle,
	// they include the vertices and the edge nodes on the RT element,
	// in counter-clockwise progression
	var ii int
	for n := 0; n < 3; n++ {
		EdgeX[ii], EdgeY[ii] = dfr.VX.DataP[n], dfr.VY.DataP[n]
		ii++
		offset := 2*NpInt + n*NpEdge
		for i := 0; i < NpEdge; i++ {
			EdgeX[ii], EdgeY[ii] = dfr.FluxX.DataP[offset+i], dfr.FluxY.DataP[offset+i]
			ii++
		}
	}
	gm = geometry2D.TriangulateTriangle(EdgeX, EdgeY,
		dfr.FluxX.DataP[:NpInt], dfr.FluxY.DataP[:NpInt])
	return
}

func TestEdgeProjectionIsentropicVortex(t *testing.T) {
	iv := isentropic_vortex.NewIVortex(5, 0, 0, 1.4, 0)
	NMin := 2
	NMax := 2
	for N := NMin; N <= NMax; N++ {
		angle := 140.
		dfr := CreateEquiTriMesh(N, angle)
		var (
			Np = dfr.FluxElement.Np
		)
		ep := dfr.FluxElement.projectInteriorToEdges()
		// Q := utils.NewMatrix(Np, 4)
		for i := 0; i < Np; i++ {
			x, y := dfr.FluxX.DataP[i], dfr.FluxY.DataP[i]
			rho, rhoU, rhoV, rhoE := iv.GetStateC(0, x, y)
			fmt.Printf("Q[%.2f,%.2f] = [%.2f,%.2f,%.2f,%.2f]\n",
				x, y, rho, rhoU, rhoV, rhoE)
		}
		// Let's use edge 2 to test the edge interpolator
		ei := NewEdgeInterpolator(ep.s[2], 2., 5.5)
		vMin, vMax := 1., 1.1691
		samples := []float64{vMax, vMax, vMax, vMin} // Mach 1.2 density
		p := ei.FitAndBoundPolynomial(samples, vMin, vMax)
		for i, c := range p.Coeffs {
			fmt.Printf("Coeff[%d]:%.2f\n", i, c)
		}
		for i, s := range ep.s[2] {
			fmt.Printf("%d: [s,Sample]:[%5.4f,%5.4f]:p:%5.4f\n",
				i, s, samples[i], p.Evaluate(s))
		}

		samples = []float64{vMax, vMin, vMax, vMax} // Mach 1.2 density
		p = ei.FitAndBoundPolynomial(samples, vMin, vMax)
		for i, c := range p.Coeffs {
			fmt.Printf("Coeff[%d]:%.2f\n", i, c)
		}
		for i, s := range ep.s[2] {
			fmt.Printf("%d: [s,Sample]:[%5.4f,%5.4f]:p:%5.4f\n",
				i, s, samples[i], p.Evaluate(s))
		}
	}
}

func TestEdgeProjectionShockCapture(t *testing.T) {
	NMin := 2
	NMax := 2
	for N := NMin; N <= NMax; N++ {
		angle := 140.
		dfr := CreateEquiTriMesh(N, angle)
		ep := dfr.FluxElement.projectInteriorToEdges()
		fmt.Println(ep.indices[2])
		fmt.Println(ep.s[2])
		// Let's use edge 2 to test the edge interpolator
		ei := NewEdgeInterpolator(ep.s[2], 2., 5.5)
		vMin, vMax := 1., 1.1691
		samples := []float64{vMax, vMax, vMax, vMin} // Mach 1.2 density
		p := ei.FitAndBoundPolynomial(samples, vMin, vMax)
		for i, c := range p.Coeffs {
			fmt.Printf("Coeff[%d]:%.2f\n", i, c)
		}
		for i, s := range ep.s[2] {
			fmt.Printf("%d: [s,Sample]:[%5.4f,%5.4f]:p:%5.4f\n",
				i, s, samples[i], p.Evaluate(s))
		}

		samples = []float64{vMax, vMin, vMax, vMax} // Mach 1.2 density
		p = ei.FitAndBoundPolynomial(samples, vMin, vMax)
		for i, c := range p.Coeffs {
			fmt.Printf("Coeff[%d]:%.2f\n", i, c)
		}
		for i, s := range ep.s[2] {
			fmt.Printf("%d: [s,Sample]:[%5.4f,%5.4f]:p:%5.4f\n",
				i, s, samples[i], p.Evaluate(s))
		}
	}
}

type EdgeProjection struct {
	indices [3][]int
	r, s    [3][]float64
}

// Sort sorts the three fields together for each of the three groups.
// The sorting key is chosen as follows:
//   - For group 0 and group 1, sort by r.
//   - For group 2, sort by s.
func (ep *EdgeProjection) Sort() {
	// Loop over each of the 3 groups.
	for group := 0; group < 3; group++ {
		n := len(ep.indices[group])
		// Build a temporary slice of entries holding the three values.
		type entry struct {
			idx int
			r   float64
			s   float64
		}
		tmp := make([]entry, n)
		for i := 0; i < n; i++ {
			tmp[i] = entry{
				idx: ep.indices[group][i],
				r:   ep.r[group][i],
				s:   ep.s[group][i],
			}
		}
		// Choose the key based on group:
		// Group 2 sorts by s; groups 0 and 1 sort by r.
		if group == 2 {
			sort.Slice(tmp, func(i, j int) bool {
				return tmp[i].s < tmp[j].s
			})
		} else {
			sort.Slice(tmp, func(i, j int) bool {
				return tmp[i].r < tmp[j].r
			})
		}
		// Copy sorted entries back into the three fields.
		for i := 0; i < n; i++ {
			ep.indices[group][i] = tmp[i].idx
			ep.r[group][i] = tmp[i].r
			ep.s[group][i] = tmp[i].s
		}
	}
}

func (rt *RTElement) projectInteriorToEdges() (ep EdgeProjection) {
	ep = EdgeProjection{}
	NpInt := rt.NpInt
	NpEdge := rt.NpEdge
	tangentNeg := func(v [2]float64) (tt [2]float64) {
		// Left hand rule tangent to follow counter-clockwise edge traversal
		tt = [2]float64{-v[1], v[0]}
		tNorm := math.Sqrt(tt[0]*tt[0] + tt[1]*tt[1])
		tt[0] /= tNorm
		tt[1] /= tNorm
		return
	}
	rsFromTangentProj := func(r, s float64,
		normal, origin [2]float64) (rP, sP float64) {
		tt := tangentNeg(normal)
		dot := tt[0]*(r-origin[0]) + tt[1]*(s-origin[1])
		rP = tt[0]*dot + origin[0]
		sP = tt[1]*dot + origin[1]
		return
	}
	fmt.Printf("Number of edge points: %d, Interior Points: %d\n",
		NpEdge, NpInt)
	var origin [2]float64
	// Indices used for projection of interior points, filtered from total
	for nEdge := 0; nEdge < 3; nEdge++ {
		switch nEdge {
		case 0:
			origin = [2]float64{-1, -1}
		case 1:
			origin = [2]float64{1, -1}
		case 2:
			origin = [2]float64{-1, 1}
		}
		fmt.Printf("nEdge: %d\n", nEdge)
		offset := 2*NpInt + nEdge*NpEdge
		// Compose a set of interior points, then filter them to find the
		// optimal group for projection
		var points [][2]float64
		for i := 0; i < NpInt; i++ {
			r, s := rt.R.AtVec(i), rt.S.AtVec(i)
			// fmt.Printf("point%d = [%5.4f,%5.4f]\n", i, r, s)
			points = append(points, [2]float64{r, s})
		}
		// for i := 2 * NpInt; i < 2*NpInt+3*NpEdge; i++ {
		// 	r, s := rt.R.AtVec(i), rt.S.AtVec(i)
		// fmt.Printf("edge point%d = [%5.4f,%5.4f]\n", i, r, s)
		// }
		edgeNormal := rt.DOFVectors[offset] // Edge normal vector
		ep.indices[nEdge] = filterInteriorPointIndices(points,
			edgeNormal.Eval(), origin)
		// fmt.Println("indices: ", indices[nEdge])
		var rP, sP float64
		for _, i := range ep.indices[nEdge] {
			rr, ss := rt.R.AtVec(i), rt.S.AtVec(i)
			b_j := rt.DOFVectors[offset] // Edge normal vector
			rP, sP = rsFromTangentProj(rr, ss, b_j.Eval(), origin)
			ep.r[nEdge] = append(ep.r[nEdge], rP)
			ep.s[nEdge] = append(ep.s[nEdge], sP)
		}
		ep.Sort()
		for ii, i := range ep.indices[nEdge] {
			rr, ss := rt.R.AtVec(i), rt.S.AtVec(i)
			rP = ep.r[nEdge][ii]
			sP = ep.s[nEdge][ii]
			fmt.Printf("IntPt:%d Proj(%5.4f,%5.4f)=[%5.4f,%5.4f]\n",
				i, rr, ss, rP, sP)
		}
	}
	return
}

func filterInteriorPointIndices(points [][2]float64, normal, origin [2]float64) []int {
	// Tolerance for considering projected points as colocated.
	tol := 0.01

	// Pre-compute the tangent (using the left-hand rule) once.
	tangent := [2]float64{-normal[1], normal[0]}
	tNorm := math.Sqrt(tangent[0]*tangent[0] + tangent[1]*tangent[1])
	tangent[0] /= tNorm
	tangent[1] /= tNorm

	// rsFromTangentProj projects a point (r,s) along the tangent direction.
	rsFromTangentProj := func(r, s float64) (rP, sP float64) {
		// dot product with the tangent (shifted by origin)
		dot := tangent[0]*(r-origin[0]) + tangent[1]*(s-origin[1])
		rP = tangent[0]*dot + origin[0]
		sP = tangent[1]*dot + origin[1]
		return
	}

	// Compute the perpendicular (normal) distance of a point from the edge.
	normalDistance := func(r, s float64) float64 {
		return (r-origin[0])*normal[0] + (s-origin[1])*normal[1]
	}

	// Start with every point indexed.
	indices := make([]int, len(points))
	for i := range points {
		indices[i] = i
	}

	for {
		// Map from the quantized projected coordinate to the list of indices that share it.
		projMap := make(map[[2]float64][]int)
		for _, idx := range indices {
			p := points[idx]
			rP, sP := rsFromTangentProj(p[0], p[1])
			// Quantize the projection using the tolerance
			key := [2]float64{math.Round(rP/tol) * tol, math.Round(sP/tol) * tol}
			projMap[key] = append(projMap[key], idx)
		}

		// From each group, select the point that is closest to the edge (in absolute distance).
		newIndices := make([]int, 0, len(indices))
		for _, group := range projMap {
			if len(group) == 1 {
				newIndices = append(newIndices, group[0])
			} else {
				bestIdx := group[0]
				bestDist := math.Abs(normalDistance(points[bestIdx][0], points[bestIdx][1]))
				for _, idx := range group[1:] {
					d := math.Abs(normalDistance(points[idx][0], points[idx][1]))
					if d < bestDist {
						bestIdx = idx
						bestDist = d
					}
				}
				newIndices = append(newIndices, bestIdx)
			}
		}

		// Check for stability: if the set of indices is unchanged (both in count and content), we're done.
		if len(newIndices) == len(indices) {
			sort.Ints(newIndices)
			sortedOld := make([]int, len(indices))
			copy(sortedOld, indices)
			sort.Ints(sortedOld)
			equal := true
			for i := range newIndices {
				if newIndices[i] != sortedOld[i] {
					equal = false
					break
				}
			}
			if equal {
				return newIndices
			}
		}
		// Otherwise, continue filtering with the new set.
		indices = newIndices
	}
}

func TestEdgeInterpolation(t *testing.T) {
	NMin := 7
	NMax := 7
	// for N := NMin; N <= NMax; N++ {
	for N := NMin; N <= NMax; N++ {
		angle := 140.
		dfr := CreateEquiTriMesh(N, angle)
		NpInt := dfr.FluxElement.NpInt
		NpFlux := dfr.FluxElement.Np
		// NEdge := dfr.FluxElement.NpEdge
		// fmt.Printf("NInt[%d] = %d, NInt/3 = %d, Remainder=%d, NEdge=%d\n",
		// 	N, NpInt, NpInt/3, NpInt%3, NEdge)
		M := GradientMassMatrix(dfr)
		QSol := utils.NewMatrix(NpInt, 4)
		QFlux := utils.NewMatrix(NpFlux, 4)
		field := setShockConditions(dfr.FluxX.DataP, 2, 0)
		for n := 0; n < 4; n++ {
			for i := 0; i < NpFlux; i++ {
				if i < NpInt {
					QSol.Set(i, n, field[i+n*NpFlux])
				}
				QFlux.Set(i, n, field[i+n*NpFlux])
			}
		}
		// QSol.Print("QSol")
		// QFlux.Print("QFlux")
		MInv := M.Transpose().Mul(M).InverseWithCheck().Mul(M.Transpose())
		// MInv.Print("MInv")
		Dens := CopyDensityFromQAndEdges(dfr, QSol, QFlux)
		Grad := MInv.Mul(Dens)
		// Grad.Print("Grad")
		GradA, GradB, GradC := Grad.DataP[0], Grad.DataP[1], Grad.DataP[2]
		n1, n2 := GradB, GradC
		norm := math.Sqrt(n1*n1 + n2*n2)
		n1 /= norm
		n2 /= norm
		fmt.Printf("Order[%d] - Shock normal vector: [%5.2f,%5.2f]\n",
			N, n1, n2)
		// Now compute an estimate of the shock threshold rho_thresh as the
		// average of min/max density in the element
		rhoMin, rhoMax := Dens.Min(), Dens.Max()
		rho_thresh := (rhoMin + rhoMax) / 2.
		RhoGrad := func(x, y float64) (RhoFit float64) {
			// Quickly determines which side of the shock surface we are on
			RhoFit = GradA + GradB*x + GradC*y
			return
		}
		// If RhoGrad > rho_thresh, we're on the post shock side
		_, _ = rho_thresh, RhoGrad
		var leftRight [2][]int
		for i := 0; i < NpInt; i++ {
			x, y := dfr.SolutionX.At(i, 0), dfr.SolutionY.At(i, 0)
			rhoFit := RhoGrad(x, y)
			var label string
			if rhoFit > rho_thresh {
				label = "post"
				leftRight[1] = append(leftRight[1], i)
			} else {
				label = "pre"
				leftRight[0] = append(leftRight[0], i)
			}
			fmt.Printf("Solution Point[%d][%5.2f,"+
				"%5.2f] %s shock rho=[%5.2f], rhoFit=[%5.2f]\n",
				i, x, y, label, QSol.At(i, 0), rhoFit)
		}
		if testing.Verbose() {
			PlotShockPoints(dfr, leftRight, 30*time.Second)
		}
		// edgeGroups := groupInteriorPoints(dfr.FluxX.DataP, dfr.FluxY.DataP,
		// 	NpInt, NEdge)
		// fmt.Println(edgeGroups)
		// if testing.Verbose() {
		// 	PlotEdgeGroups(dfr, edgeGroups)
		// }
	}
}

func PlotShockPoints(dfr *DFR2D, leftRight [2][]int, delay time.Duration) {
	var Line []float32
	addVertex := func(x, y, width float64) {
		wo2 := float32(width / 2.)
		xl, yl := float32(x), float32(y)
		Line = append(Line, xl-wo2, yl, xl+wo2, yl)
		Line = append(Line, xl, yl-wo2, xl, yl+wo2)
	}
	addLineSegment := func(x1, y1, x2, y2 float64) {
		x1l, y1l := float32(x1), float32(y1)
		x2l, y2l := float32(x2), float32(y2)
		Line = append(Line, x1l, y1l, x2l, y2l)
	}

	ch := chart2d.NewChart2D(-1, 1, -1, 1, 1024, 1024, utils2.WHITE,
		utils2.BLACK, 0.8)
	Line = append(Line, 0, -10, 0, 10)
	ch.AddLine(Line, utils2.WHITE)
	Line = []float32{}
	// Np := dfr.FluxElement.Np
	for k := 0; k < dfr.K; k++ {
		for ii, position := range leftRight {
			for _, pt := range position {
				addVertex(dfr.SolutionX.At(pt, k), dfr.SolutionY.At(pt, k), 0.025)
			}
			switch ii {
			case 0:
				ch.AddLine(Line, utils2.RED)
			case 1:
				ch.AddLine(Line, utils2.GREEN)
			}
			Line = []float32{}
		}
		// for i := 0; i < Np; i++ {
		// 	addVertex(dfr.FluxX.At(i, k), dfr.FluxY.At(i, k), 0.025)
		// }
	}
	for k := 0; k < dfr.K; k++ {
		verts := dfr.Tris.GetTriVerts(uint32(k))
		for ii := 0; ii < 3; ii++ {
			x1, y1 := dfr.VX.DataP[verts[ii]], dfr.VY.DataP[verts[ii]]
			ip := ii + 1
			if ii == 2 {
				ip = 0
			}
			x2, y2 := dfr.VX.DataP[verts[ip]], dfr.VY.DataP[verts[ip]]
			addLineSegment(x1, y1, x2, y2)
			ch.AddLine(Line, utils2.WHITE)
			Line = []float32{}
		}
	}
	// ch.NewWindow("Unit Triangle", 0.9, screen.AUTO)
	// Line = []float32{}
	// addLineSegment(-1, -1, 1, -1)
	// addLineSegment(1, -1, -1, 1)
	// addLineSegment(-1, 1, -1, -1)
	// ch.AddLine(Line, utils2.WHITE)
	time.Sleep(delay)
}

func CopyDensityFromQAndEdges(dfr *DFR2D, QSol,
	QFlux utils.Matrix) (Dens utils.Matrix) {
	var (
		// Np     = dfr.SolutionElement.Np
		// NpFlux = dfr.FluxElement.Np
		NpInt  = dfr.FluxElement.NpInt
		NpEdge = dfr.FluxElement.NpEdge
		K      = 1
	)
	Dens = utils.NewMatrix(NpInt+3*NpEdge, K)
	for i := 0; i < NpInt; i++ {
		Dens.Set(i, 0, QSol.At(i, 0))
	}
	for i := NpInt; i < NpInt+3*NpEdge; i++ {
		Dens.Set(i, 0, QFlux.At(i+NpInt, 0))
	}
	return
}

func GradientMassMatrix(dfr *DFR2D) (M utils.Matrix) {
	var (
		// NpFlux = dfr.FluxElement.Np
		NpInt  = dfr.FluxElement.NpInt
		NpEdge = dfr.FluxElement.NpEdge
	)
	// Compute gradient using all points in flux element with density
	// First compose a mass matrix with all of the point locations
	M = utils.NewMatrix(NpInt+3*NpEdge, 3) // Points x (1:X:Y)
	for i := 0; i < NpInt; i++ {
		M.Set(i, 0, 1.)
		M.Set(i, 1, dfr.FluxX.DataP[i])
		M.Set(i, 2, dfr.FluxY.DataP[i])
	}
	for i := NpInt; i < NpInt+3*NpEdge; i++ {
		M.Set(i, 0, 1.)
		M.Set(i, 1, dfr.FluxX.DataP[i+NpInt])
		M.Set(i, 2, dfr.FluxY.DataP[i+NpInt])
	}
	// M.Print("M")
	return
}

// projectToEdge projects a 2D point (x, y) onto the edge defined by (x0, y0) to (x1, y1)
// using the edge’s tangent direction.
func projectToEdge(x, y, x0, y0, x1, y1 float64) float64 {
	// Compute the unit tangent vector of the edge.
	tanX := x1 - x0
	tanY := y1 - y0
	length := math.Hypot(tanX, tanY)
	tanX /= length
	tanY /= length

	// The 1D coordinate is the dot product of (x-x0, y-y0) with the tangent.
	return (x-x0)*tanX + (y-y0)*tanY
}

func PlotEdgeGroups(dfr *DFR2D, edgeGroups [3][]int) {
	var Line []float32
	addVertex := func(x, y, width float64) {
		wo2 := float32(width / 2.)
		xl, yl := float32(x), float32(y)
		Line = append(Line, xl-wo2, yl, xl+wo2, yl)
		Line = append(Line, xl, yl-wo2, xl, yl+wo2)
	}
	addLineSegment := func(x1, y1, x2, y2 float64) {
		x1l, y1l := float32(x1), float32(y1)
		x2l, y2l := float32(x2), float32(y2)
		Line = append(Line, x1l, y1l, x2l, y2l)
	}

	ch := chart2d.NewChart2D(-1, 1, -1, 1, 1024, 1024, utils2.WHITE,
		utils2.BLACK, 0.8)
	// Np := dfr.FluxElement.Np
	for k := 0; k < dfr.K; k++ {
		for ii, edge := range edgeGroups {
			for _, pt := range edge {
				addVertex(dfr.FluxX.At(pt, k), dfr.FluxY.At(pt, k), 0.025)
			}
			switch ii {
			case 0:
				ch.AddLine(Line, utils2.RED)
			case 1:
				ch.AddLine(Line, utils2.GREEN)
			case 2:
				ch.AddLine(Line, utils2.BLUE)
			}
			Line = []float32{}
		}
		// for i := 0; i < Np; i++ {
		// 	addVertex(dfr.FluxX.At(i, k), dfr.FluxY.At(i, k), 0.025)
		// }
	}
	for k := 0; k < dfr.K; k++ {
		verts := dfr.Tris.GetTriVerts(uint32(k))
		for ii := 0; ii < 3; ii++ {
			x1, y1 := dfr.VX.DataP[verts[ii]], dfr.VY.DataP[verts[ii]]
			ip := ii + 1
			if ii == 2 {
				ip = 0
			}
			x2, y2 := dfr.VX.DataP[verts[ip]], dfr.VY.DataP[verts[ip]]
			addLineSegment(x1, y1, x2, y2)
			switch ii {
			case 0:
				ch.AddLine(Line, utils2.RED)
			case 1:
				ch.AddLine(Line, utils2.GREEN)
			case 2:
				ch.AddLine(Line, utils2.BLUE)
			}
			Line = []float32{}
		}
	}
	// ch.NewWindow("Unit Triangle", 0.9, screen.AUTO)
	// Line = []float32{}
	// addLineSegment(-1, -1, 1, -1)
	// addLineSegment(1, -1, -1, 1)
	// addLineSegment(-1, 1, -1, -1)
	// ch.AddLine(Line, utils2.WHITE)
	time.Sleep(30 * time.Second)
}

// candidate holds per‐interior point distance info.
type candidate struct {
	index     int
	distances [3]float64 // min distance to each edge group.
	primary   int        // index of closest edge.
	secondary int        // index of next closest edge.
	delta     float64    // difference between secondary and primary distance.
}

// groupInteriorPoints first assigns each interior point (indices 0..NpInt-1)
// to the edge with the minimum distance, then rebalances if any group exceeds its target.
// The edges are at:
//
//	Edge 1: indices [2*NpInt, 2*NpInt+NpEdge)
//	Edge 2: indices [2*NpInt+NpEdge, 2*NpInt+2*NpEdge)
//	Edge 3: indices [2*NpInt+2*NpEdge, 2*NpInt+3*NpEdge)
func groupInteriorPoints(x, y []float64, NpInt, NpEdge int) [3][]int {
	cands := make([]candidate, NpInt)
	// squaredDistance returns the squared Euclidean distance.
	squaredDistance := func(x1, y1, x2, y2 float64) float64 {
		dx := x1 - x2
		dy := y1 - y2
		return dx*dx + dy*dy
	}
	// For each interior point, compute its distance to each edge.
	for i := 0; i < NpInt; i++ {
		best := math.MaxFloat64
		second := math.MaxFloat64
		bestEdge := -1
		secondEdge := -1

		// For each of the three edges.
		for e := 0; e < 3; e++ {
			start := 2*NpInt + e*NpEdge
			end := start + NpEdge
			distE := math.MaxFloat64
			for j := start; j < end; j++ {
				d := squaredDistance(x[i], y[i], x[j], y[j])
				if d < distE {
					distE = d
				}
			}
			cands[i].distances[e] = distE
			// Determine best and second best edges.
			if distE < best {
				second = best
				secondEdge = bestEdge
				best = distE
				bestEdge = e
			} else if distE < second {
				second = distE
				secondEdge = e
			}
		}
		cands[i].index = i
		cands[i].primary = bestEdge
		cands[i].secondary = secondEdge
		cands[i].delta = second - best
	}

	// Initial assignment: each interior point goes to its primary edge.
	assigned := make([]int, NpInt)
	for i := 0; i < NpInt; i++ {
		assigned[i] = cands[i].primary
	}

	// Count assignments per edge.
	counts := [3]int{}
	for i := 0; i < NpInt; i++ {
		counts[assigned[i]]++
	}

	// Determine target counts.
	base := NpInt / 3
	rem := NpInt % 3
	targets := [3]int{}
	for e := 0; e < 3; e++ {
		if e < rem {
			targets[e] = base + 1
		} else {
			targets[e] = base
		}
	}

	// Rebalancing: for any point in an over-target group, if its secondary candidate is under-target,
	// reassign it. We do a simple greedy pass.
	for {
		moved := false
		for i := 0; i < NpInt; i++ {
			curGroup := assigned[i]
			// Only consider points in groups that have too many points.
			if counts[curGroup] > targets[curGroup] {
				sec := cands[i].secondary
				// If the secondary candidate exists and is under target, reassign.
				if sec >= 0 && counts[sec] < targets[sec] {
					assigned[i] = sec
					counts[curGroup]--
					counts[sec]++
					moved = true
					// Break out to restart the loop after a move.
					break
				}
			}
		}
		if !moved {
			break
		}
	}

	// Build final groups.
	var groups [3][]int
	for i, grp := range assigned {
		groups[grp] = append(groups[grp], i)
	}
	return groups
}

func TestPlotEquiTri(t *testing.T) {
	// dfr := CreateEquiTriMesh(1, 0.15)
	// dfr := CreateEquiTriMesh(2, 0.15)
	dfr := CreateEquiTriMesh(2, 0.15)
	_ = dfr
	// if testing.Verbose() {
	// 	PlotDFRElements(dfr)
	// }
	Mach := 1.8
	Alpha := 10.
	fmt.Printf("Testing Mach %5.2f, Shock Angle:%5.2f\n", Mach, Alpha)
	U1, U2 := ShockConditions(Mach, Alpha)
	PrintU("U1", U1)
	PrintU("U2", U2)
	K := 1
	Np := dfr.SolutionElement.Np
	Q := utils.NewMatrix(Np, 4*K)
	for i := 0; i < Np; i++ {
		x, _ := dfr.SolutionX.At(i, 0), dfr.SolutionY.At(i, 0)
		// fmt.Printf("x,y[%d] = [%5.2f,%5.2f]\n", i, x, y)
		for n := 0; n < 4; n++ {
			if x < 0 {
				Q.DataP[n+i*4*K] = U1[n]
			} else {
				Q.DataP[n+i*4*K] = U2[n]
			}
		}
	}
	// Q.Print("Q1")
	// dfr.FluxEdgeInterp.Print("FluxEdgeInterp")
	QEdge := dfr.FluxEdgeInterp.Mul(Q)
	// QEdge.Print("QEdge")
	// dfr.FluxElement.ProjectFunctionOntoDOF()
	NpFluxEdge := dfr.FluxElement.NpEdge
	var CorrectU [4]float64
	var RMS, Max, MaxPercent [4]float64
	for i := 0; i < NpFluxEdge*3; i++ {
		x, _ := dfr.FluxX.At(i, 0), dfr.FluxY.At(i, 0)
		// fmt.Printf("x,y[%d] = [%5.2f,%5.2f]\n", i, x, y)
		if x < 0 {
			CorrectU = U1
		} else {
			CorrectU = U2
		}
		for n := 0; n < 4; n++ {
			err := QEdge.At(i, n) - CorrectU[n]
			RMS[n] += err * err
			if math.Abs(err) > Max[n] {
				Max[n] = err
				if math.Abs(CorrectU[n]) < 0.00001 {
					MaxPercent[n] = 100. * err
				} else {
					MaxPercent[n] = 100. * err / CorrectU[n]
				}
			}
		}
	}
	for n := 0; n < 4; n++ {
		RMS[n] = math.Sqrt(RMS[n] / float64(NpFluxEdge*3))
		fmt.Printf("Error[%d]:RMS%+5.2f, Max:%+5.2f,%+5.1f%%\n",
			n, RMS[n], Max[n], MaxPercent[n])
	}
}

func PrintU(label string, U [4]float64) {
	fmt.Printf("%s:", label)
	for i := 0; i < 4; i++ {
		fmt.Printf("%0.2f ", U[i])
	}
	fmt.Printf("\n")
}

func TestMachConditions(t *testing.T) {
	U1, U2 := ShockConditions(2, 0)
	fmt.Printf("U1: %5.2f,%5.2f,%5.2f,%5.2f\n", U1[0], U1[1], U1[2], U1[3])
	fmt.Printf("U2: %5.2f,%5.2f,%5.2f,%5.2f\n", U2[0], U2[1], U2[2], U2[3])
	U1, U2 = ShockConditions(2, 15)
	fmt.Printf("U1: %5.2f,%5.2f,%5.2f,%5.2f\n", U1[0], U1[1], U1[2], U1[3])
	fmt.Printf("U2: %5.2f,%5.2f,%5.2f,%5.2f\n", U2[0], U2[1], U2[2], U2[3])
}

// ShockConditions computes post-shock properties given pre-shock Mach number (M1) and shock angle (alpha)
// where alpha = 0 corresponds to a normal shock.
func ShockConditions(M1, alpha float64) (U1, U2 [4]float64) {
	// Alpha = angle of incoming flow relative to the shock front (0 = normal shock).
	var (
		gamma = 1.4
	)

	// Convert angles to radians
	alphaRad := alpha * math.Pi / 180.0

	// Pre-shock conditions (non-dimensional)
	rho1 := 1.0
	p1 := 1.0
	u1 := M1
	v1 := 0.0 // Flow is initially aligned along x (free stream)

	// Resolve incoming velocity into normal/tangential components relative to the shock
	u1n := u1*math.Cos(alphaRad) + v1*math.Sin(alphaRad)  // Normal velocity component
	u1t := -u1*math.Sin(alphaRad) + v1*math.Cos(alphaRad) // Tangential component (unchanged)

	// Normal Mach number
	M1n := u1n / math.Sqrt(gamma*p1/rho1)

	// If subsonic normal component, no shock forms
	if M1n < 1.0 {
		En1 := p1/(gamma-1) + 0.5*rho1*(u1*u1+v1*v1)
		U1 = [4]float64{rho1, rho1 * u1, rho1 * v1, En1}
		U2 = U1
		return
	}

	// Rankine-Hugoniot relations (normal direction)
	rhoRatio := ((gamma + 1) * M1n * M1n) / ((gamma-1)*M1n*M1n + 2)
	pRatio := 1 + (2*gamma/(gamma+1))*(M1n*M1n-1)
	u2n_to_u1n := ((gamma-1)*M1n*M1n + 2) / ((gamma + 1) * M1n * M1n)

	// Post-shock normal values
	rho2 := rho1 * rhoRatio
	p2 := p1 * pRatio
	u2n := u1n * u2n_to_u1n

	// Tangential velocity is unchanged across the shock (slip condition for inviscid flow)
	u2t := u1t

	// Now transform back to x,y coordinates (undo rotation)
	u2 := u2n*math.Cos(alphaRad) - u2t*math.Sin(alphaRad)
	v2 := u2n*math.Sin(alphaRad) + u2t*math.Cos(alphaRad)

	// Total energy (consistent with post-shock pressure and velocity)
	En1 := p1/(gamma-1) + 0.5*rho1*(u1*u1+v1*v1)
	En2 := p2/(gamma-1) + 0.5*rho2*(u2*u2+v2*v2)

	// Pack conservative variables into U1 and U2
	U1 = [4]float64{rho1, rho1 * u1, rho1 * v1, En1}
	U2 = [4]float64{rho2, rho2 * u2, rho2 * v2, En2}

	return
}

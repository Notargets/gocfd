package DG2D

import (
	"math"
	"testing"
)

func TestInterpolationMethodComparison(t *testing.T) {
	// Focus on order 3 for detailed analysis
	NMax := 5
	for N := 3; N < NMax; N++ {
		t.Logf("\n============================================================")
		t.Logf("ORDER %d INTERPOLATION METHOD COMPARISON", N)
		t.Logf("============================================================")

		// Create DG element (interior 2D, order N)
		le := NewLagrangeElement2D(N)
		// Create RT element (edges 1D, order N+1)
		rt := NewRTElement(N+1, SimplexRTBasis, OptimizedEdgePoints)

		t.Logf("DG element: %d interior points", le.Np)
		t.Logf("RT element: %d interior points, %d points per edge", rt.NpInt, rt.NpEdge)

		// Test multiple polynomial types to expose different numerical behaviors
		testCases := []struct {
			name        string
			poly        func(r, s float64) float64
			description string
		}{
			{
				name: "smooth_polynomial",
				poly: func(r, s float64) float64 {
					// Smooth polynomial using most of the available order
					result := 0.0
					for i := 0; i <= N-1; i++ { // Use N-1 to stay well within bounds
						for j := 0; j <= N-1-i; j++ {
							coeff := 1.0 + 0.1*float64(i) + 0.05*float64(j)
							result += coeff * math.Pow(r, float64(i)) * math.Pow(s, float64(j))
						}
					}
					return result
				},
				description: "Smooth polynomial using order N-1 terms",
			},
			{
				name: "oscillatory_polynomial",
				poly: func(r, s float64) float64 {
					// More oscillatory - might expose interpolation issues
					return 1.0 + math.Sin(math.Pi*r) + math.Cos(math.Pi*s) + 0.5*math.Sin(2*math.Pi*r*s)
				},
				description: "Oscillatory trigonometric function",
			},
			{
				name: "sharp_gradient",
				poly: func(r, s float64) float64 {
					// Sharp but smooth gradient - mimics shock-like behavior
					sigma := 0.3
					x0, y0 := 0.2, 0.2
					return 1.0 + math.Tanh((math.Sqrt((r-x0)*(r-x0)+(s-y0)*(s-y0))-0.5)/sigma)
				},
				description: "Sharp gradient (tanh profile)",
			},
			{
				name: "boundary_layer",
				poly: func(r, s float64) float64 {
					// Boundary layer-like profile
					return 1.0 + math.Exp(-10*(r+1)*(r+1)) + 0.5*math.Exp(-8*(s+1)*(s+1))
				},
				description: "Boundary layer-like exponential profile",
			},
		}

		// Test each polynomial type
		for _, testCase := range testCases {
			t.Logf("\n--- Testing %s: %s ---", testCase.name, testCase.description)

			// Evaluate test polynomial at DG interior points
			interiorValues := make([]float64, le.Np)
			for i := 0; i < le.Np; i++ {
				r, s := le.R.AtVec(i), le.S.AtVec(i)
				interiorValues[i] = testCase.poly(r, s)
			}

			minVal, maxVal := minFloat64(interiorValues), maxFloat64(interiorValues)
			dynamicRange := maxVal - minVal
			t.Logf("Interior values range: [%.6f, %.6f], dynamic range: %.6f",
				minVal, maxVal, dynamicRange)

			// Test each edge
			edgeResults := make([]EdgeMethodComparison, 3)
			for edgeNum := 0; edgeNum < 3; edgeNum++ {
				edgeResults[edgeNum] = compareInterpolationMethods(
					le, rt, edgeNum, interiorValues, testCase.poly, t)
			}

			// Summarize results for this polynomial
			summarizePolynomialResults(testCase.name, edgeResults, t)
		}
	}
}

type EdgeMethodComparison struct {
	EdgeNum            int
	DirectError        float64
	DirectConditioning float64
	ModalError         float64
	ModalConditioning  float64
	DirectSuccess      bool
	ModalSuccess       bool
}

func compareInterpolationMethods(le *LagrangeElement2D, rt *RTElement, edgeNum int,
	interiorValues []float64, testPoly func(float64, float64) float64,
	t *testing.T) EdgeMethodComparison {

	t.Logf("\n  Edge %d:", edgeNum)

	// Method 1: Direct interpolation
	directValues, directCond, directErr := directInterpolationMethod(
		le, rt, edgeNum, interiorValues, testPoly)

	// Method 2: Modal transfer (using the WORKING implementation)
	modalValues, modalCond, modalErr := modalTransferMethodFixed(
		le, rt, edgeNum, interiorValues, testPoly)

	// Compare the actual interpolated values between methods
	methodDifference := 0.0
	if len(directValues) == len(modalValues) {
		for i := 0; i < len(directValues); i++ {
			diff := math.Abs(directValues[i] - modalValues[i])
			if diff > methodDifference {
				methodDifference = diff
			}
		}
	}

	// Define success criteria (more lenient for higher orders)
	errorTol := 1e-8 // Relaxed from machine precision
	condTol := 1e16  // Allow high conditioning but not extreme

	directSuccess := directErr < errorTol && directCond < condTol
	modalSuccess := modalErr < errorTol && modalCond < condTol

	t.Logf("    Direct: error=%.2e, cond=%.2e, success=%v",
		directErr, directCond, directSuccess)
	t.Logf("    Modal:  error=%.2e, cond=%.2e, success=%v",
		modalErr, modalCond, modalSuccess)
	t.Logf("    Method difference: %.2e (max |direct - modal|)", methodDifference)

	// Compare methods
	if modalSuccess && !directSuccess {
		t.Logf("    ✓ Modal transfer WINS (works where direct fails)")
	} else if directSuccess && !modalSuccess {
		t.Logf("    ⚠ Direct interpolation WINS (modal fails)")
	} else if modalSuccess && directSuccess {
		if modalCond < directCond/10 {
			t.Logf("    ✓ Modal transfer WINS (much better conditioning)")
		} else if directCond < modalCond/10 {
			t.Logf("    ✓ Direct interpolation WINS (much better conditioning)")
		} else {
			t.Logf("    = Both methods work (similar performance)")
		}

		// Also show if the methods give very different results
		if methodDifference > errorTol {
			t.Logf("    ⚠ WARNING: Methods disagree significantly (diff=%.2e)", methodDifference)
		}
	} else {
		t.Logf("    ✗ Both methods fail")
	}

	return EdgeMethodComparison{
		EdgeNum:            edgeNum,
		DirectError:        directErr,
		DirectConditioning: directCond,
		ModalError:         modalErr,
		ModalConditioning:  modalCond,
		DirectSuccess:      directSuccess,
		ModalSuccess:       modalSuccess,
	}
}

// func directInterpolationMethod(le *LagrangeElement2D, rt *RTElement, edgeNum int,
// 	interiorValues []float64, testPoly func(float64, float64) float64) (
// 	[]float64, float64, float64) {
//
// 	// Get edge points
// 	edgePoints := getEdgePoints(rt, edgeNum)
//
// 	// Build interpolation matrix
// 	RFlux := make([]float64, len(edgePoints))
// 	SFlux := make([]float64, len(edgePoints))
// 	for i, pt := range edgePoints {
// 		RFlux[i], SFlux[i] = pt.R, pt.S
// 	}
//
// 	RFluxVec := utils.NewVector(len(RFlux), RFlux)
// 	SFluxVec := utils.NewVector(len(SFlux), SFlux)
//
// 	Vedge := le.JB2D.Vandermonde2D(le.N, RFluxVec, SFluxVec)
// 	InterpMatrix := Vedge.Mul(le.JB2D.Vinv)
//
// 	// Apply interpolation
// 	interiorMatrix := utils.NewMatrix(len(interiorValues), 1, interiorValues)
// 	edgeMatrix := InterpMatrix.Mul(interiorMatrix)
//
// 	edgeValues := make([]float64, len(edgePoints))
// 	for i := 0; i < len(edgePoints); i++ {
// 		edgeValues[i] = edgeMatrix.At(i, 0)
// 	}
//
// 	// Compute conditioning and error
// 	conditioning := InterpMatrix.ConditionNumber()
// 	maxError := computeMaxError(edgePoints, edgeValues, testPoly)
//
// 	return edgeValues, conditioning, maxError
// }
//
// func modalTransferMethodFixed(le *LagrangeElement2D, rt *RTElement, edgeNum int,
// 	interiorValues []float64, testPoly func(float64, float64) float64) (
// 	[]float64, float64, float64) {
//
// 	// Step 1: Convert to modal coefficients
// 	interiorMatrix := utils.NewMatrix(len(interiorValues), 1, interiorValues)
// 	modalCoeffs2D := le.JB2D.Vinv.Mul(interiorMatrix)
//
// 	// Step 2: Build modal transfer matrix using the WORKING implementation
// 	transferMatrix := buildActualModalTransferMatrix(le, rt, edgeNum)
//
// 	// Step 3: Transfer to 1D edge modes
// 	modalCoeffs1D := transferMatrix.Mul(modalCoeffs2D)
//
// 	// Step 4: Convert to nodal values on edge using your Jacobi polynomial functions
// 	edgePoints := getEdgePoints(rt, edgeNum)
// 	edge1DCoords := make([]float64, len(edgePoints))
// 	for i, pt := range edgePoints {
// 		edge1DCoords[i] = mapToEdge1D(pt.R, pt.S, edgeNum)
// 	}
//
// 	// Build 1D Vandermonde using your JacobiBasis1D
// 	edgeOrder := rt.NpEdge - 1
// 	edge1DCoordVec := utils.NewVector(len(edge1DCoords), edge1DCoords)
//
// 	// Create 1D Jacobi basis (using Legendre polynomials: alpha=0, beta=0)
// 	jb1d := NewJacobiBasis1D(edgeOrder, edge1DCoordVec, 0.0, 0.0)
// 	Vedge1D := jb1d.Vandermonde1D()
//
// 	edgeMatrix := Vedge1D.Mul(modalCoeffs1D)
//
// 	edgeValues := make([]float64, len(edgePoints))
// 	for i := 0; i < len(edgePoints); i++ {
// 		edgeValues[i] = edgeMatrix.At(i, 0)
// 	}
//
// 	// Compute conditioning and error
// 	conditioning := transferMatrix.ConditionNumber()
// 	maxError := computeMaxError(edgePoints, edgeValues, testPoly)
//
// 	return edgeValues, conditioning, maxError
// }
//
// func buildActualModalTransferMatrix(le *LagrangeElement2D, rt *RTElement, edgeNum int) utils.Matrix {
// 	num2DModes := le.Np
// 	num1DModes := rt.NpEdge
//
// 	transferMatrix := utils.NewMatrix(num1DModes, num2DModes)
//
// 	// For each 2D mode, compute its projection onto each 1D mode by numerical integration
// 	for j := 0; j < num2DModes; j++ {
// 		order2D := le.JB2D.Order2DAtJ[j]
// 		i, k := order2D[0], order2D[1]
//
// 		for m := 0; m < num1DModes; m++ {
// 			// Compute projection by numerical integration along the edge
// 			coeff := computeActualModalTransferCoeff(le, i, k, m, edgeNum)
// 			transferMatrix.Set(m, j, coeff)
// 		}
// 	}
//
// 	return transferMatrix
// }
//
// func computeActualModalTransferCoeff(le *LagrangeElement2D, i, k int, m int, edgeNum int) float64 {
// 	// Compute projection by numerical integration: ∫[edge] Phi_2D(r,s) * Phi_1D(xi) dxi
// 	// Need to account for normalization differences between 2D and 1D bases
//
// 	// Use high-order Gauss-Legendre quadrature
// 	nQuad := 15
// 	quadPoints, quadWeights := getGaussLegendreQuadrature(nQuad)
//
// 	integral := 0.0
//
// 	for q := 0; q < nQuad; q++ {
// 		xi := quadPoints[q]
// 		weight := quadWeights[q]
//
// 		// Map xi to (r,s) on the edge
// 		r, s := mapEdgeParamToRS(xi, edgeNum)
//
// 		// Evaluate 2D basis function Phi_{i,k}(r,s) using YOUR actual implementation
// 		RVec := utils.NewVector(1, []float64{r})
// 		SVec := utils.NewVector(1, []float64{s})
// 		basis2D := le.JB2D.Simplex2DP(RVec, SVec, i, k)
// 		val2D := basis2D[0]
//
// 		// Evaluate 1D Legendre polynomial P_m(xi)
// 		val1D := evaluateLegendrePolynomial(xi, m)
//
// 		// Add to integral (Jacobian = 1 for our parametrization)
// 		integral += weight * val2D * val1D
// 	}
//
// 	// Apply normalization correction
// 	// Your 2D basis includes √2 factor, 1D Legendre has √(2/(2m+1)) normalization
// 	normalization1D := math.Sqrt(2.0 / (2.0*float64(m) + 1.0))
//
// 	return integral / normalization1D
// }
//
// func getGaussLegendreQuadrature(n int) ([]float64, []float64) {
// 	// High-precision Gauss-Legendre quadrature points and weights
// 	switch n {
// 	case 15:
// 		points := []float64{
// 			-0.9879925180204854, -0.9372733924007059, -0.8482065834104272, -0.7244177313601700,
// 			-0.5709721726085388, -0.3941513470775634, -0.2011940939974345, 0.0000000000000000,
// 			0.2011940939974345, 0.3941513470775634, 0.5709721726085388, 0.7244177313601700,
// 			0.8482065834104272, 0.9372733924007059, 0.9879925180204854,
// 		}
// 		weights := []float64{
// 			0.0307532419961173, 0.0703660474881081, 0.1071592204671719, 0.1395706779261543,
// 			0.1662692058169939, 0.1861610000155622, 0.1984314853271116, 0.2025782419255613,
// 			0.1984314853271116, 0.1861610000155622, 0.1662692058169939, 0.1395706779261543,
// 			0.1071592204671719, 0.0703660474881081, 0.0307532419961173,
// 		}
// 		return points, weights
// 	default:
// 		// Fallback to 5-point rule
// 		points := []float64{-0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831, 0.9061798459386640}
// 		weights := []float64{0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891}
// 		return points, weights
// 	}
// }
//
// func mapEdgeParamToRS(xi float64, edgeNum int) (float64, float64) {
// 	switch edgeNum {
// 	case 0: // Bottom edge: s = -1, r = xi
// 		return xi, -1.0
// 	case 1: // Hypotenuse: r + s = 0, r = xi, s = -xi
// 		return xi, -xi
// 	case 2: // Left edge: r = -1, s = xi
// 		return -1.0, xi
// 	default:
// 		return 0.0, 0.0
// 	}
// }
//
// func evaluateLegendrePolynomial(x float64, n int) float64 {
// 	if n == 0 {
// 		return 1.0
// 	} else if n == 1 {
// 		return x
// 	}
//
// 	p0, p1 := 1.0, x
// 	for k := 2; k <= n; k++ {
// 		p2 := (float64(2*k-1)*x*p1 - float64(k-1)*p0) / float64(k)
// 		p0, p1 = p1, p2
// 	}
// 	return p1
// }

func summarizePolynomialResults(polyName string, results []EdgeMethodComparison, t *testing.T) {
	t.Logf("\n  Summary for %s:", polyName)

	directWins := 0
	modalWins := 0
	bothWork := 0
	bothFail := 0

	for _, result := range results {
		if result.DirectSuccess && result.ModalSuccess {
			bothWork++
		} else if result.DirectSuccess && !result.ModalSuccess {
			directWins++
		} else if !result.DirectSuccess && result.ModalSuccess {
			modalWins++
		} else {
			bothFail++
		}
	}

	t.Logf("    Direct only: %d edges, Modal only: %d edges", directWins, modalWins)
	t.Logf("    Both work: %d edges, Both fail: %d edges", bothWork, bothFail)

	if modalWins > directWins {
		t.Logf("    → Modal transfer shows better robustness")
	} else if directWins > modalWins {
		t.Logf("    → Direct interpolation shows better robustness")
	} else {
		t.Logf("    → Methods show similar robustness")
	}
}

// Helper functions
// func computeMaxError(edgePoints []EdgePoint, computed []float64,
// 	exact func(float64, float64) float64) float64 {
// 	maxError := 0.0
// 	for i, pt := range edgePoints {
// 		exactVal := exact(pt.R, pt.S)
// 		error := math.Abs(computed[i] - exactVal)
// 		if error > maxError {
// 			maxError = error
// 		}
// 	}
// 	return maxError
// }
//
// func getEdgePoints(rt *RTElement, edgeNum int) []EdgePoint {
// 	NpInt := rt.NpInt
// 	points := make([]EdgePoint, rt.NpEdge)
// 	startIdx := 2*NpInt + edgeNum*rt.NpEdge
//
// 	for i := 0; i < rt.NpEdge; i++ {
// 		idx := startIdx + i
// 		points[i] = EdgePoint{
// 			R: rt.R.AtVec(idx),
// 			S: rt.S.AtVec(idx),
// 		}
// 	}
// 	return points
// }
//
// func mapToEdge1D(r, s float64, edgeNum int) float64 {
// 	switch edgeNum {
// 	case 0: // Bottom edge
// 		return r
// 	case 1: // Hypotenuse
// 		return r
// 	case 2: // Left edge
// 		return s
// 	default:
// 		return 0.0
// 	}
// }

func minFloat64(arr []float64) float64 {
	min := arr[0]
	for _, v := range arr {
		if v < min {
			min = v
		}
	}
	return min
}

func maxFloat64(arr []float64) float64 {
	max := arr[0]
	for _, v := range arr {
		if v > max {
			max = v
		}
	}
	return max
}

// type EdgePoint struct {
// 	R, S float64
// }

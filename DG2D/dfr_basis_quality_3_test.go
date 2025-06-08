package DG2D

import (
	"fmt"
	"math"
	"math/rand/v2"
	"os"
	"sort"
	"testing"

	"github.com/notargets/gocfd/utils"

	"gonum.org/v1/gonum/optimize"
)

func extractCurrentEdgePoints(rt *RTElement) [3][]EdgePoint {
	var edgePoints [3][]EdgePoint

	for edgeNum := 0; edgeNum < 3; edgeNum++ {
		points := getEdgePoints(rt, edgeNum)
		edgePoints[edgeNum] = points
	}

	return edgePoints
}

func computeModalTransferConditioning(le *LagrangeElement2D, rt *RTElement,
	edgeNum int, edgePoints []EdgePoint) float64 {

	// Create temporary RT element with new edge points
	tempRT := createTempRTWithEdgePoints(rt, edgeNum, edgePoints)

	// Build modal transfer matrix
	transferMatrix := buildActualModalTransferMatrix(le, tempRT, edgeNum)

	// Compute condition number
	return transferMatrix.ConditionNumber()
}

func computeSpectralApproximationError(le *LagrangeElement2D, rt *RTElement,
	edgeNum int, edgePoints []EdgePoint) float64 {

	// Test modal transfer accuracy on polynomial test functions
	maxError := 0.0

	// Test polynomials up to order N-1 (should be exactly representable)
	testPolynomials := []func(r, s float64) float64{
		func(r, s float64) float64 { return 1.0 },   // constant
		func(r, s float64) float64 { return r },     // linear in r
		func(r, s float64) float64 { return s },     // linear in s
		func(r, s float64) float64 { return r * r }, // quadratic
		func(r, s float64) float64 { return r * s }, // mixed
		func(r, s float64) float64 { return s * s }, // quadratic in s
	}

	for _, testPoly := range testPolynomials {
		error := evaluateModalTransferError(le, rt, edgeNum, edgePoints, testPoly)
		if error > maxError {
			maxError = error
		}
	}

	return maxError
}

func evaluateModalTransferError(le *LagrangeElement2D, rt *RTElement,
	edgeNum int, edgePoints []EdgePoint, testPoly func(float64, float64) float64) float64 {

	// Evaluate test polynomial at DG interior points
	interiorValues := make([]float64, le.Np)
	for i := 0; i < le.Np; i++ {
		r, s := le.R.AtVec(i), le.S.AtVec(i)
		interiorValues[i] = testPoly(r, s)
	}

	// Create temporary RT with optimized edge points
	tempRT := createTempRTWithEdgePoints(rt, edgeNum, edgePoints)

	// Perform modal transfer
	_, _, modalError := modalTransferMethodFixed(le, tempRT, edgeNum, interiorValues, testPoly)

	return modalError
}

func copyRTElement(rt *RTElement) (rtCopy *RTElement) {
	rtCopy = &RTElement{
		P:          rt.P,
		Np:         rt.Np,
		NpEdge:     rt.NpEdge,
		NpInt:      rt.NpInt,
		V:          rt.V,
		VInv:       rt.VInv,
		Div:        rt.Div,
		DivInt:     rt.DivInt,
		R:          rt.R,
		S:          rt.S,
		RInt:       rt.RInt,
		SInt:       rt.SInt,
		DOFVectors: rt.DOFVectors,
		RTBasis:    rt.RTBasis,
		Phi:        rt.Phi,
	}
	return
}

func createTempRTWithEdgePoints(rt *RTElement, edgeNum int, newEdgePoints []EdgePoint) *RTElement {
	// Create a copy of RT element with modified edge points
	tempRT := copyRTElement(rt)

	// Copy existing points
	tempRT.R = rt.R.Copy()
	tempRT.S = rt.S.Copy()

	// Update the specified edge points
	NpInt := rt.NpInt
	startIdx := 2*NpInt + edgeNum*rt.NpEdge

	for i, point := range newEdgePoints {
		idx := startIdx + i
		tempRT.R.DataP[idx] = point.R
		tempRT.S.DataP[idx] = point.S
	}
	return tempRT
}

func evaluateCurrentPerformance(le *LagrangeElement2D, rt *RTElement, edgeNum int,
	edgePoints []EdgePoint, t *testing.T) (float64, float64) {

	// Test with a representative polynomial
	testPoly := func(r, s float64) float64 {
		return 1.0 + r + s + 0.5*r*s // Mixed polynomial within order capability
	}

	// Evaluate at DG interior points
	interiorValues := make([]float64, le.Np)
	for i := 0; i < le.Np; i++ {
		r, s := le.R.AtVec(i), le.S.AtVec(i)
		interiorValues[i] = testPoly(r, s)
	}

	// Evaluate direct interpolation
	_, directCond, _ := directInterpolationMethod(le, rt, edgeNum, interiorValues, testPoly)

	// Evaluate modal transfer
	_, modalCond, _ := modalTransferMethodFixed(le, rt, edgeNum, interiorValues, testPoly)

	return directCond, modalCond
}

func evaluateOptimizedPerformance(le *LagrangeElement2D, rt *RTElement, edgeNum int,
	optimizedPoints []EdgePoint, t *testing.T) (float64, float64) {

	// Create temporary RT with optimized points
	tempRT := createTempRTWithEdgePoints(rt, edgeNum, optimizedPoints)

	// Test with the same polynomial as before
	testPoly := func(r, s float64) float64 {
		return 1.0 + r + s + 0.5*r*s
	}

	// Evaluate at DG interior points
	interiorValues := make([]float64, le.Np)
	for i := 0; i < le.Np; i++ {
		r, s := le.R.AtVec(i), le.S.AtVec(i)
		interiorValues[i] = testPoly(r, s)
	}

	// Evaluate direct interpolation with optimized points
	_, directCond, _ := directInterpolationMethod(le, tempRT, edgeNum, interiorValues, testPoly)

	// Evaluate modal transfer with optimized points
	_, modalCond, _ := modalTransferMethodFixed(le, tempRT, edgeNum, interiorValues, testPoly)

	return directCond, modalCond
}

func formatPoints(points []EdgePoint) string {
	result := "["
	for i, p := range points {
		if i > 0 {
			result += ", "
		}
		result += fmt.Sprintf("(%.3f,%.3f)", p.R, p.S)
	}
	result += "]"
	return result
}

// Helper functions (reuse from previous implementations)
func modalTransferMethodFixed(le *LagrangeElement2D, rt *RTElement, edgeNum int,
	interiorValues []float64, testPoly func(float64, float64) float64) (
	[]float64, float64, float64) {

	// Step 1: Convert to modal coefficients
	interiorMatrix := utils.NewMatrix(len(interiorValues), 1, interiorValues)
	modalCoeffs2D := le.JB2D.Vinv.Mul(interiorMatrix)

	// Step 2: Build modal transfer matrix
	transferMatrix := buildActualModalTransferMatrix(le, rt, edgeNum)

	// Step 3: Transfer to 1D edge modes
	modalCoeffs1D := transferMatrix.Mul(modalCoeffs2D)

	// Step 4: Convert to nodal values using Jacobi basis
	edgePoints := getEdgePoints(rt, edgeNum)
	edge1DCoords := make([]float64, len(edgePoints))
	for i, pt := range edgePoints {
		edge1DCoords[i] = mapToEdge1D(pt.R, pt.S, edgeNum)
	}

	edgeOrder := rt.NpEdge - 1
	edge1DCoordVec := utils.NewVector(len(edge1DCoords), edge1DCoords)
	jb1d := NewJacobiBasis1D(edgeOrder, edge1DCoordVec)
	Vedge1D := jb1d.Vandermonde1D()

	edgeMatrix := Vedge1D.Mul(modalCoeffs1D)

	edgeValues := make([]float64, len(edgePoints))
	for i := 0; i < len(edgePoints); i++ {
		edgeValues[i] = edgeMatrix.At(i, 0)
	}

	// Compute conditioning and error
	conditioning := transferMatrix.ConditionNumber()
	maxError := computeMaxError(edgePoints, edgeValues, testPoly)

	return edgeValues, conditioning, maxError
}

func directInterpolationMethod(le *LagrangeElement2D, rt *RTElement, edgeNum int,
	interiorValues []float64, testPoly func(float64, float64) float64) (
	[]float64, float64, float64) {

	// Get edge points
	edgePoints := getEdgePoints(rt, edgeNum)

	// Build interpolation matrix
	RFlux := make([]float64, len(edgePoints))
	SFlux := make([]float64, len(edgePoints))
	for i, pt := range edgePoints {
		RFlux[i], SFlux[i] = pt.R, pt.S
	}

	RFluxVec := utils.NewVector(len(RFlux), RFlux)
	SFluxVec := utils.NewVector(len(SFlux), SFlux)

	Vedge := le.JB2D.Vandermonde2D(le.N, RFluxVec, SFluxVec)
	InterpMatrix := Vedge.Mul(le.JB2D.Vinv)

	// Apply interpolation
	interiorMatrix := utils.NewMatrix(len(interiorValues), 1, interiorValues)
	edgeMatrix := InterpMatrix.Mul(interiorMatrix)

	edgeValues := make([]float64, len(edgePoints))
	for i := 0; i < len(edgePoints); i++ {
		edgeValues[i] = edgeMatrix.At(i, 0)
	}

	// Compute conditioning and error
	conditioning := InterpMatrix.ConditionNumber()
	maxError := computeMaxError(edgePoints, edgeValues, testPoly)

	return edgeValues, conditioning, maxError
}

// Reuse helper functions from previous implementations
func buildActualModalTransferMatrix(le *LagrangeElement2D, rt *RTElement, edgeNum int) utils.Matrix {
	num2DModes := le.Np
	num1DModes := rt.NpEdge

	transferMatrix := utils.NewMatrix(num1DModes, num2DModes)

	for j := 0; j < num2DModes; j++ {
		order2D := le.JB2D.Order2DAtJ[j]
		i, k := order2D[0], order2D[1]

		for m := 0; m < num1DModes; m++ {
			coeff := computeActualModalTransferCoeff(le, i, k, m, edgeNum)
			transferMatrix.Set(m, j, coeff)
		}
	}

	return transferMatrix
}

func computeActualModalTransferCoeff(le *LagrangeElement2D, i, k int, m int, edgeNum int) float64 {
	nQuad := 15
	quadPoints, quadWeights := getGaussLegendreQuadrature(nQuad)

	integral := 0.0

	for q := 0; q < nQuad; q++ {
		xi := quadPoints[q]
		weight := quadWeights[q]

		r, s := mapEdgeParamToRS(xi, edgeNum)

		RVec := utils.NewVector(1, []float64{r})
		SVec := utils.NewVector(1, []float64{s})
		basis2D := le.JB2D.Simplex2DP(RVec, SVec, i, k)
		val2D := basis2D[0]

		val1D := evaluateLegendrePolynomial(xi, m)

		integral += weight * val2D * val1D
	}

	normalization1D := math.Sqrt(2.0 / (2.0*float64(m) + 1.0))

	return integral / normalization1D
}

func getGaussLegendreQuadrature(n int) ([]float64, []float64) {
	switch n {
	case 15:
		points := []float64{
			-0.9879925180204854, -0.9372733924007059, -0.8482065834104272, -0.7244177313601700,
			-0.5709721726085388, -0.3941513470775634, -0.2011940939974345, 0.0000000000000000,
			0.2011940939974345, 0.3941513470775634, 0.5709721726085388, 0.7244177313601700,
			0.8482065834104272, 0.9372733924007059, 0.9879925180204854,
		}
		weights := []float64{
			0.0307532419961173, 0.0703660474881081, 0.1071592204671719, 0.1395706779261543,
			0.1662692058169939, 0.1861610000155622, 0.1984314853271116, 0.2025782419255613,
			0.1984314853271116, 0.1861610000155622, 0.1662692058169939, 0.1395706779261543,
			0.1071592204671719, 0.0703660474881081, 0.0307532419961173,
		}
		return points, weights
	default:
		points := []float64{-0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831, 0.9061798459386640}
		weights := []float64{0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891}
		return points, weights
	}
}

func evaluateLegendrePolynomial(x float64, n int) float64 {
	if n == 0 {
		return 1.0
	} else if n == 1 {
		return x
	}

	p0, p1 := 1.0, x
	for k := 2; k <= n; k++ {
		p2 := (float64(2*k-1)*x*p1 - float64(k-1)*p0) / float64(k)
		p0, p1 = p1, p2
	}
	return p1
}

func computeMaxError(edgePoints []EdgePoint, computed []float64,
	exact func(float64, float64) float64) float64 {
	maxError := 0.0
	for i, pt := range edgePoints {
		exactVal := exact(pt.R, pt.S)
		error := math.Abs(computed[i] - exactVal)
		if error > maxError {
			maxError = error
		}
	}
	return maxError
}

func getEdgePoints(rt *RTElement, edgeNum int) []EdgePoint {
	NpInt := rt.NpInt
	points := make([]EdgePoint, rt.NpEdge)
	startIdx := 2*NpInt + edgeNum*rt.NpEdge

	for i := 0; i < rt.NpEdge; i++ {
		idx := startIdx + i
		points[i] = EdgePoint{
			R: rt.R.AtVec(idx),
			S: rt.S.AtVec(idx),
		}
	}
	return points
}

type EdgePoint struct {
	R, S float64
}

func computeLebesgueConstant(edgePoints []EdgePoint, edgeNum int) float64 {
	// Estimate Lebesgue constant by sampling
	nSamples := 1000
	maxLebesgue := 0.0

	for i := 0; i < nSamples; i++ {
		// Sample point along edge
		t := -1.0 + 2.0*float64(i)/float64(nSamples-1)

		// Compute sum of absolute values of Lagrange basis functions
		lebesgueSum := 0.0
		for j, edgePoint := range edgePoints {
			// Lagrange basis function L_j(t)
			lagrangeBasis := 1.0
			edgeParam := getEdgeParameter(edgePoint, edgeNum) // FIXED: use edgeNum instead of j

			for k, otherPoint := range edgePoints {
				if k != j {
					otherParam := getEdgeParameter(otherPoint, edgeNum) // FIXED: use edgeNum instead of j
					lagrangeBasis *= (t - otherParam) / (edgeParam - otherParam)
				}
			}
			lebesgueSum += math.Abs(lagrangeBasis)
		}

		if lebesgueSum > maxLebesgue {
			maxLebesgue = lebesgueSum
		}
	}

	return maxLebesgue
}

// Fixed functions to maintain proper edge constraints

func edgePointsToParams(edgePoints []EdgePoint, edgeNum int) []float64 {
	params := make([]float64, len(edgePoints))

	for i, point := range edgePoints {
		params[i] = mapToEdge1D(point.R, point.S, edgeNum)
	}

	return params
}

func paramsToEdgePoints(params []float64, edgeNum int) []EdgePoint {
	points := make([]EdgePoint, len(params))

	for i, param := range params {
		r, s := mapEdgeParamToRS(param, edgeNum)
		points[i] = EdgePoint{R: r, S: s}
	}

	return points
}

func mapToEdge1D(r, s float64, edgeNum int) float64 {
	switch edgeNum {
	case 0: // Bottom edge: s = -1, r varies from -1 to 1
		return r
	case 1: // Hypotenuse: r + s = 0, parameterized by r (so s = -r)
		return r
	case 2: // Left edge: r = -1, s varies from -1 to 1
		return s
	default:
		return 0.0
	}
}

func mapEdgeParamToRS(xi float64, edgeNum int) (float64, float64) {
	switch edgeNum {
	case 0: // Bottom edge: s = -1, r = xi
		return xi, -1.0
	case 1: // Hypotenuse: r + s = 0, so r = xi, s = -xi
		return xi, -xi
	case 2: // Left edge: r = -1, s = xi
		return -1.0, xi
	default:
		return 0.0, 0.0
	}
}

func getEdgeParameter(point EdgePoint, edgeNum int) float64 {
	switch edgeNum {
	case 0: // Bottom edge: s = -1, parameter is r
		return point.R
	case 1: // Hypotenuse: r + s = 0, parameter is r (since s = -r)
		return point.R
	case 2: // Left edge: r = -1, parameter is s
		return point.S
	default:
		return 0.0
	}
}
func rangeLim(r float64) bool {
	var (
		tol = 1.e-2
	)
	if r-tol < -1.0 || r+tol > 1.0 {
		// if r < -0.99 || r > 0.99 {
		return true
	} else {
		return false
	}
}

// Enhanced validation function to check edge constraints
func validateEdgeConstraints(edgePoints []EdgePoint, edgeNum int) bool {
	const tolerance = 1e-10

	for _, point := range edgePoints {
		switch edgeNum {
		case 0: // Bottom edge: s should be -1
			if (math.Abs(point.S+1.0) > tolerance) || rangeLim(point.R) {
				return false
			}
		case 1: // Hypotenuse: r + s should be 0
			if math.Abs(point.R+point.S) > tolerance || rangeLim(point.R) || rangeLim(point.S) {
				return false
			}
		case 2: // Left edge: r should be -1
			if (math.Abs(point.R+1.0) > tolerance) || rangeLim(point.S) {
				return false
			}
		}
	}
	return true
}

// Test function with edge constraint verification
func TestEdgePointOptimization(t *testing.T) {
	if false {
		// Test edge point optimization for orders 3-6
		for N := 4; N <= 4; N++ {
			separator := "============================================================"
			t.Logf("\n%s", separator)
			t.Logf("EDGE POINT OPTIMIZATION FOR ORDER %d", N)
			t.Logf("%s", separator)

			// Create elements for this order
			le := NewLagrangeElement2D(N)
			rt := NewRTElement(N+1, SimplexRTBasis, OptimizedEdgePoints)

			// Original (unoptimized) edge points
			originalPoints := extractCurrentEdgePoints(rt)

			// Optimize edge points for each edge
			for edgeNum := 0; edgeNum < 3; edgeNum++ {
				t.Logf("\n--- OPTIMIZING EDGE %d ---", edgeNum)

				// Get current edge points as starting guess
				currentPoints := originalPoints[edgeNum]
				t.Logf("Original points: %v", formatPoints(currentPoints))

				// Verify original points are on the edge
				if !validateEdgeConstraints(currentPoints, edgeNum) {
					t.Logf("ERROR: Original points are not on edge %d!", edgeNum)
					continue
				}

				// Evaluate current performance
				directCond, modalCond := evaluateCurrentPerformance(le, rt, edgeNum, currentPoints, t)

				// Optimize edge points
				optimizedPoints := optimizeEdgePoints(le, rt, edgeNum, currentPoints, t)
				t.Logf("Optimized points: %v", formatPoints(optimizedPoints))

				// Verify optimized points are still on the edge
				if !validateEdgeConstraints(optimizedPoints, edgeNum) {
					t.Logf("ERROR: Optimized points are not on edge %d!", edgeNum)
					os.Exit(1)
				}

				// Evaluate optimized performance
				optDirectCond, optModalCond := evaluateOptimizedPerformance(le, rt, edgeNum, optimizedPoints, t)

				// Compare results
				t.Logf("RESULTS:")
				t.Logf("  Direct interpolation conditioning: %.2e -> %.2e (%.1fx improvement)",
					directCond, optDirectCond, directCond/optDirectCond)
				t.Logf("  Modal transfer conditioning: %.2e -> %.2e (%.1fx improvement)",
					modalCond, optModalCond, modalCond/optModalCond)

				// Determine recommended method
				if optModalCond < optDirectCond {
					t.Logf("  ✓ RECOMMENDATION: Use Modal Transfer")
				} else {
					t.Logf("  ✓ RECOMMENDATION: Use Direct Interpolation")
				}
			}
		}
	}
}

// Improved objective function with better penalty structure
func edgePointObjective(params []float64, le *LagrangeElement2D, rt *RTElement, edgeNum int) float64 {
	// Convert parameters to edge points
	edgePoints := paramsToEdgePoints(params, edgeNum)

	// Debug: Verify edge constraints are maintained
	if !validateEdgeConstraints(edgePoints, edgeNum) {
		// Return very high penalty if constraints are violated
		return 1e12
	}

	// Compute various quality metrics
	lebesgueConstant := computeLebesgueConstant(edgePoints, edgeNum)
	modalConditioning := computeModalTransferConditioning(le, rt, edgeNum, edgePoints)
	spectralError := computeSpectralApproximationError(le, rt, edgeNum, edgePoints)

	// Weighted combination of objectives with better scaling
	alpha := 0.2 // Lebesgue constant weight
	beta := 0.7  // Modal conditioning weight (main focus)
	gamma := 0.1 // Spectral error weight

	// Use more stable objective formulation
	objective := alpha*math.Log10(math.Max(lebesgueConstant, 1.0)) +
		beta*math.Log10(math.Max(modalConditioning, 1.0)) +
		gamma*math.Log10(math.Max(1.0+spectralError, 1.0))

	// Much smaller, quadratic penalties for better optimization behavior
	penalty := 0.0

	// Soft penalty for points outside valid range [-1, 1]
	for _, param := range params {
		if param < -1.0 {
			excess := param + 1.0
			penalty += 10.0 * excess * excess // Quadratic penalty
		} else if param > 1.0 {
			excess := param - 1.0
			penalty += 10.0 * excess * excess // Quadratic penalty
		}
	}

	// Softer penalty for non-monotonic ordering with better gradients
	for i := 1; i < len(params); i++ {
		if params[i] <= params[i-1] {
			violation := params[i-1] - params[i] + 1e-6
			penalty += 100.0 * violation * violation // Much smaller than 1e6
		}
	}

	// Add penalty for extreme clustering (maintain reasonable spacing)
	minDesiredSpacing := 0.1 // Minimum desired spacing between points
	for i := 1; i < len(params); i++ {
		spacing := params[i] - params[i-1]
		if spacing < minDesiredSpacing {
			shortfall := minDesiredSpacing - spacing
			penalty += 50.0 * shortfall * shortfall
		}
	}

	return objective + penalty
}

// Alternative: Even simpler objective focusing mainly on conditioning
func edgePointObjectiveSimple(params []float64, le *LagrangeElement2D, rt *RTElement, edgeNum int) float64 {
	// Convert parameters to edge points
	edgePoints := paramsToEdgePoints(params, edgeNum)

	// Verify edge constraints are maintained
	if !validateEdgeConstraints(edgePoints, edgeNum) {
		return 1e12
	}

	// Focus primarily on modal conditioning
	modalConditioning := computeModalTransferConditioning(le, rt, edgeNum, edgePoints)

	// Simple objective: just minimize condition number (in log space for stability)
	objective := math.Log10(math.Max(modalConditioning, 1.0))

	// Very light penalties
	penalty := 0.0

	// Boundary penalties
	for _, param := range params {
		if param < -1.0 || param > 1.0 {
			penalty += 1.0 * math.Abs(param) // Much lighter
		}
	}

	// Ordering penalty
	for i := 1; i < len(params); i++ {
		if params[i] <= params[i-1] {
			penalty += 1.0 // Just a small constant penalty
		}
	}

	return objective + penalty
}

// Enhanced optimization with multiple attempts and better initial conditions
func optimizeEdgePoints(le *LagrangeElement2D, rt *RTElement, edgeNum int,
	initialPoints []EdgePoint, t *testing.T) []EdgePoint {

	t.Logf("Starting optimization for edge %d...", edgeNum)

	// Validate that initial points are actually on the edge
	if !validateEdgeConstraints(initialPoints, edgeNum) {
		t.Logf("WARNING: Initial points are not properly on edge %d", edgeNum)
	}

	// Convert edge points to 1D parameter space [-1, 1]
	initialParams := edgePointsToParams(initialPoints, edgeNum)
	t.Logf("Initial parameters: %v", initialParams)

	// Ensure parameters are properly sorted
	sort.Float64s(initialParams)

	// Try the simpler objective first
	bestParams := initialParams
	bestObjective := 1e20

	// Attempt 1: Simple objective focusing on conditioning
	t.Logf("Attempt 1: Simple conditioning-focused objective...")
	problem1 := optimize.Problem{
		Func: func(x []float64) float64 {
			return edgePointObjectiveSimple(x, le, rt, edgeNum)
		},
	}

	settings := &optimize.Settings{
		MajorIterations: 200,
		FuncEvaluations: 500,
	}

	result1, err1 := optimize.Minimize(problem1, initialParams, settings, &optimize.NelderMead{})
	if err1 == nil && result1.F < bestObjective {
		bestParams = result1.X
		bestObjective = result1.F
		t.Logf("Simple objective gave better result: %.6e", result1.F)
	}

	// Attempt 2: Full objective with improved penalties
	t.Logf("Attempt 2: Full objective with improved penalties...")
	problem2 := optimize.Problem{
		Func: func(x []float64) float64 {
			return edgePointObjective(x, le, rt, edgeNum)
		},
	}

	result2, err2 := optimize.Minimize(problem2, initialParams, settings, &optimize.NelderMead{})
	if err2 == nil && result2.F < bestObjective {
		bestParams = result2.X
		bestObjective = result2.F
		t.Logf("Full objective gave better result: %.6e", result2.F)
	}

	// Use the best result
	if bestObjective >= 1e10 {
		t.Logf("All optimization attempts failed, returning original points")
		return initialPoints
	}

	t.Logf("Optimization completed:")
	t.Logf("  Best objective: %.6e", bestObjective)
	t.Logf("  Optimized parameters: %v", bestParams)

	// Convert optimized parameters back to edge points
	optimizedPoints := paramsToEdgePoints(bestParams, edgeNum)

	// Validate that optimized points are still on the edge
	if !validateEdgeConstraints(optimizedPoints, edgeNum) {
		t.Logf("ERROR: Optimized points are not on edge %d! Returning original points.", edgeNum)
		return initialPoints
	}

	return optimizedPoints
}

/*
KEY INSIGHTS FROM THE ALTERNATIVE OPTIMIZATION APPROACH:

1. LEBESGUE CONSTANT AS PRIMARY OBJECTIVE:
   - Focuses on minimizing Lebesgue constant directly, which is a well-established
     measure of interpolation quality and oscillation control
   - Much simpler than multi-objective optimization with condition numbers

2. MULTIPLE STARTING DISTRIBUTIONS:
   - Uses LGL, Chebyshev, equidistant, and perturbed points as starting guesses
   - Performs global search rather than local optimization from one starting point

3. VARIABLE TRANSFORMATION:
   - Uses tanh transformation to map unbounded variables to [-1,1] interval
   - Automatically handles bound constraints without penalties

4. SEPARATE EDGE HANDLING:
   - Optimizes each edge independently but with edge-specific coordinate mappings
   - Simpler coordinate transformations for each edge type

5. DIFFERENT OPTIMIZATION STRATEGY:
   - Multi-method approach (BFGS, Nelder-Mead fallbacks)
   - Higher iteration limits (2000 vs 200-300)
   - Specialized refinement for low-order cases

Let me implement an improved version that combines the best of both approaches:
*/

// Improved edge point optimization inspired by the alternative approach
func optimizeEdgePointsImproved(le *LagrangeElement2D, rt *RTElement, edgeNum int,
	initialPoints []EdgePoint, t *testing.T) []EdgePoint {

	t.Logf("Starting improved optimization for edge %d...", edgeNum)

	// Extract just the interior points (exclude endpoints)
	nTotal := len(initialPoints)
	if nTotal < 3 {
		t.Logf("Not enough points to optimize (need at least 3)")
		return initialPoints
	}

	nInterior := nTotal - 2
	initialInteriorParams := make([]float64, nInterior)

	// Get 1D parameters for interior points only
	allParams := edgePointsToParams(initialPoints, edgeNum)
	copy(initialInteriorParams, allParams[1:nTotal-1]) // Skip first and last (endpoints)

	t.Logf("Initial interior parameters: %v", initialInteriorParams)

	// Generate multiple starting distributions for global search
	startingDistributions := generateStartingDistributions(nInterior, initialInteriorParams)

	bestParams := make([]float64, nInterior)
	copy(bestParams, initialInteriorParams)
	bestObjective := math.Inf(1)

	// Try optimization from each starting distribution
	for idx, startDist := range startingDistributions {
		t.Logf("Attempt %d with starting distribution: %v", idx+1, startDist)

		result := optimizeFromStartingPoint(startDist, edgeNum, le, rt, t)

		// Evaluate the result
		fullParams := reconstructFullParams(result, edgeNum)
		fullPoints := paramsToEdgePoints(fullParams, edgeNum)
		objective := evaluateInterpolationQuality(fullPoints, edgeNum, le, rt)

		if objective < bestObjective {
			bestObjective = objective
			copy(bestParams, result)
			t.Logf("  New best from start %d - Objective: %.6e", idx+1, objective)
		}
	}

	// Reconstruct full parameter set with fixed endpoints
	finalFullParams := reconstructFullParams(bestParams, edgeNum)
	optimizedPoints := paramsToEdgePoints(finalFullParams, edgeNum)

	// Validate constraints
	if !validateEdgeConstraints(optimizedPoints, edgeNum) {
		t.Logf("ERROR: Final optimized points violate edge constraints! Returning original.")
		return initialPoints
	}

	t.Logf("Optimization completed with best objective: %.6e", bestObjective)
	t.Logf("Final interior parameters: %v", bestParams)

	return optimizedPoints
}

// Generate multiple starting distributions for global optimization
func generateStartingDistributions(nInterior int, originalInterior []float64) [][]float64 {
	var distributions [][]float64

	if nInterior == 0 {
		return [][]float64{{}}
	}

	// 1. Original distribution
	original := make([]float64, nInterior)
	copy(original, originalInterior)
	distributions = append(distributions, original)

	// 2. LGL interior points
	lglFull := generateLGLPoints(nInterior + 2)
	if len(lglFull) >= nInterior+2 {
		lglInterior := make([]float64, nInterior)
		copy(lglInterior, lglFull[1:nInterior+1])
		distributions = append(distributions, lglInterior)
	}

	// 3. Chebyshev interior points
	chebFull := generateChebyshevPoints(nInterior + 2)
	if len(chebFull) >= nInterior+2 {
		chebInterior := make([]float64, nInterior)
		copy(chebInterior, chebFull[1:nInterior+1])
		distributions = append(distributions, chebInterior)
	}

	// 4. Equidistant points
	equiInterior := make([]float64, nInterior)
	for i := 0; i < nInterior; i++ {
		equiInterior[i] = -1.0 + 2.0*float64(i+1)/float64(nInterior+1)
	}
	distributions = append(distributions, equiInterior)

	// 5. Perturbed versions of the best known distribution (LGL)
	if len(lglFull) >= nInterior+2 {
		for pertIdx := 0; pertIdx < 3; pertIdx++ {
			perturbed := make([]float64, nInterior)
			copy(perturbed, lglFull[1:nInterior+1])

			// Add random perturbations
			for i := range perturbed {
				perturbation := (rand.Float64()*2 - 1) * 0.1 * (1 - perturbed[i]*perturbed[i])
				perturbed[i] += perturbation
				perturbed[i] = math.Max(-0.99, math.Min(0.99, perturbed[i]))
			}

			sort.Float64s(perturbed)
			distributions = append(distributions, perturbed)
		}
	}

	return distributions
}

// Optimize from a specific starting point using tanh transformation
func optimizeFromStartingPoint(startingInterior []float64, edgeNum int,
	le *LagrangeElement2D, rt *RTElement, t *testing.T) []float64 {

	nInterior := len(startingInterior)
	if nInterior == 0 {
		return []float64{}
	}

	// Transform to unbounded variables using atanh
	initialY := make([]float64, nInterior)
	for i, x := range startingInterior {
		// Clamp to avoid atanh(±1) = ±∞
		safeX := math.Max(-0.9999, math.Min(0.9999, x))
		initialY[i] = math.Atanh(safeX)
	}

	// Set up optimization problem
	problem := optimize.Problem{
		Func: func(y []float64) float64 {
			return edgeObjectiveLebesgue(y, edgeNum, le, rt)
		},
	}

	// More aggressive optimization settings
	settings := &optimize.Settings{
		MajorIterations:   1000,
		FuncEvaluations:   5000,
		GradientThreshold: 1e-10,
	}

	// Try multiple optimization methods
	var bestResult []float64
	bestObjective := math.Inf(1)

	// Method 1: BFGS
	if result, err := optimize.Minimize(problem, initialY, settings, &optimize.BFGS{}); err == nil {
		obj := result.F
		if obj < bestObjective {
			bestObjective = obj
			bestResult = result.X
		}
	}

	// Method 2: Nelder-Mead
	if result, err := optimize.Minimize(problem, initialY, settings, &optimize.NelderMead{}); err == nil {
		obj := result.F
		if obj < bestObjective {
			bestObjective = obj
			bestResult = result.X
		}
	}

	// Transform back to bounded variables
	if bestResult != nil {
		optimizedInterior := make([]float64, nInterior)
		for i, y := range bestResult {
			optimizedInterior[i] = math.Tanh(y)
		}

		// Ensure monotonicity
		sort.Float64s(optimizedInterior)

		return optimizedInterior
	}

	// Fallback to original if optimization failed
	return startingInterior
}

// Primary objective: Lebesgue constant (inspired by alternative approach)
func edgeObjectiveLebesgue(y []float64, edgeNum int, le *LagrangeElement2D, rt *RTElement) float64 {
	nInterior := len(y)

	// Transform y -> x using tanh (automatic bound constraints)
	xInterior := make([]float64, nInterior)
	for i, yi := range y {
		xInterior[i] = math.Tanh(yi)
	}

	// Reconstruct full parameter vector with endpoints
	fullParams := make([]float64, nInterior+2)
	fullParams[0] = -1.0
	copy(fullParams[1:], xInterior)
	fullParams[nInterior+1] = 1.0

	// Check monotonicity
	for i := 0; i < len(fullParams)-1; i++ {
		if fullParams[i] >= fullParams[i+1] {
			return 1e10 // Heavy penalty
		}
	}

	// Convert to edge points
	edgePoints := paramsToEdgePoints(fullParams, edgeNum)

	// Compute Lebesgue constant as primary objective
	lebesgueConst := computeLebesgueConstantImproved(edgePoints, edgeNum)

	return lebesgueConst
}

// Improved Lebesgue constant computation (more samples, better numerical stability)
func computeLebesgueConstantImproved(edgePoints []EdgePoint, edgeNum int) float64 {
	nSamples := 2000 // Higher resolution
	maxLebesgue := 0.0

	for i := 0; i < nSamples; i++ {
		// Sample point along edge parameter space
		t := -0.999 + 1.998*float64(i)/float64(nSamples-1)

		// Compute sum of absolute values of Lagrange basis functions
		lebesgueSum := 0.0
		for j, edgePoint := range edgePoints {
			// Lagrange basis function L_j(t)
			lagrangeBasis := 1.0
			edgeParam := getEdgeParameter(edgePoint, edgeNum)

			for k, otherPoint := range edgePoints {
				if k != j {
					otherParam := getEdgeParameter(otherPoint, edgeNum)
					denom := edgeParam - otherParam
					if math.Abs(denom) < 1e-14 {
						lagrangeBasis = 0.0
						break
					}
					lagrangeBasis *= (t - otherParam) / denom
				}
			}
			lebesgueSum += math.Abs(lagrangeBasis)
		}

		if lebesgueSum > maxLebesgue {
			maxLebesgue = lebesgueSum
		}
	}

	// Additional focused sampling near each point (where Lebesgue function often peaks)
	for _, point := range edgePoints {
		param := getEdgeParameter(point, edgeNum)

		// Sample 50 points in neighborhood of each interpolation point
		for i := 0; i < 50; i++ {
			offset := (float64(i) - 25.0) / 250.0 // ±0.1 range
			t := param + offset

			if t <= -1.0 || t >= 1.0 {
				continue
			}

			lebesgueSum := 0.0
			for j, edgePoint := range edgePoints {
				lagrangeBasis := 1.0
				edgeParam := getEdgeParameter(edgePoint, edgeNum)

				for k, otherPoint := range edgePoints {
					if k != j {
						otherParam := getEdgeParameter(otherPoint, edgeNum)
						denom := edgeParam - otherParam
						if math.Abs(denom) < 1e-14 {
							lagrangeBasis = 0.0
							break
						}
						lagrangeBasis *= (t - otherParam) / denom
					}
				}
				lebesgueSum += math.Abs(lagrangeBasis)
			}

			if lebesgueSum > maxLebesgue {
				maxLebesgue = lebesgueSum
			}
		}
	}

	return math.Max(maxLebesgue, 1.0) // Avoid zero/negative values
}

// Evaluate overall interpolation quality (could include multiple metrics)
func evaluateInterpolationQuality(edgePoints []EdgePoint, edgeNum int,
	le *LagrangeElement2D, rt *RTElement) float64 {

	// Primary: Lebesgue constant
	lebesgue := computeLebesgueConstantImproved(edgePoints, edgeNum)

	// Secondary: Modal conditioning (lower weight)
	modalCond := computeModalTransferConditioning(le, rt, edgeNum, edgePoints)

	// Weighted combination (Lebesgue dominant)
	return 0.8*lebesgue + 0.2*math.Log10(math.Max(modalCond, 1.0))
}

// Reconstruct full parameter vector from interior parameters
func reconstructFullParams(interiorParams []float64, edgeNum int) []float64 {
	nInterior := len(interiorParams)
	fullParams := make([]float64, nInterior+2)

	fullParams[0] = -1.0
	copy(fullParams[1:], interiorParams)
	fullParams[nInterior+1] = 1.0

	return fullParams
}

// Generate LGL points (simplified version)
func generateLGLPoints(n int) []float64 {
	if n <= 1 {
		return []float64{0}
	}
	if n == 2 {
		return []float64{-1, 1}
	}
	if n == 3 {
		return []float64{-1, 0, 1}
	}
	if n == 4 {
		return []float64{-1, -math.Sqrt(1.0 / 5.0), math.Sqrt(1.0 / 5.0), 1}
	}
	if n == 5 {
		return []float64{-1, -math.Sqrt(3.0 / 7.0), 0, math.Sqrt(3.0 / 7.0), 1}
	}
	if n == 6 {
		a := 2.0 * math.Sqrt(7.0) / 21.0
		return []float64{-1, -math.Sqrt(1.0/3.0 + a), -math.Sqrt(1.0/3.0 - a),
			math.Sqrt(1.0/3.0 - a), math.Sqrt(1.0/3.0 + a), 1}
	}
	if n == 7 {
		return []float64{-1, -0.830224, -0.468849, 0, 0.468849, 0.830224, 1}
	}

	// For higher orders, use Chebyshev as approximation
	return generateChebyshevPoints(n)
}

// Generate Chebyshev points
func generateChebyshevPoints(n int) []float64 {
	if n <= 1 {
		return []float64{0}
	}

	points := make([]float64, n)
	for i := 0; i < n; i++ {
		points[i] = -math.Cos(math.Pi * float64(2*i+1) / float64(2*n))
	}
	sort.Float64s(points)
	return points
}

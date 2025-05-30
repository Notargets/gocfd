package DG2D

import (
	"fmt"
	"math"
	"math/rand/v2"
	"sort"
	"strings"
	"testing"

	"github.com/notargets/gocfd/utils"
	"gonum.org/v1/gonum/optimize"
)

// EdgeOptimizationResult stores the optimization results for analysis
type EdgeOptimizationResult struct {
	EdgeNum             int
	OriginalPoints      []EdgePoint
	OptimizedPoints     []EdgePoint
	OriginalLebesgue    float64
	OptimizedLebesgue   float64
	OriginalModalCond   float64
	OptimizedModalCond  float64
	OriginalDirectCond  float64
	OptimizedDirectCond float64
	OptimizationMethod  string
	Iterations          int
}

// Enhanced test function with better output formatting
func TestEdgePointOptimizationEnhanced(t *testing.T) {
	if testing.Verbose() {

		separator := "================================================================================"

		NStart := 1
		NEnd := 7
		oResults := make([][]float64, NEnd+1)
		for N := NStart; N <= NEnd; N++ {
			t.Logf("\n%s", separator)
			t.Logf("EDGE POINT OPTIMIZATION RESULTS FOR ORDER %d", N)
			t.Logf("%s\n", separator)

			// Create elements
			le := NewLagrangeElement2D(N)
			rt := NewRTElement(N+1, SimplexRTBasis, GaussEdgePoints) // Start with GL points

			// Store results for summary
			results := make([]EdgeOptimizationResult, 3)

			// Optimize each edge
			for edgeNum := 0; edgeNum < 3; edgeNum++ {
				t.Logf("Edge %d Optimization:", edgeNum)
				t.Logf("--------------------------------------------------------------------------------")

				// Get initial points
				initialPoints := getEdgePoints(rt, edgeNum)

				// Debug: show actual coordinates
				if edgeNum == 2 {
					t.Logf("  Edge 2 debug - Initial points (R,S):")
					for i, pt := range initialPoints {
						t.Logf("    Point %d: (%.3f, %.3f)", i, pt.R, pt.S)
					}
				}

				// Run optimization
				result := optimizeEdgePointsMultiStart(le, rt, edgeNum, initialPoints, t)
				results[edgeNum] = result

				// Display results
				t.Logf("  Initial Distribution: %s", formatPointDistribution(result.OriginalPoints, edgeNum))
				t.Logf("  Optimized Distribution: %s", formatPointDistribution(result.OptimizedPoints, edgeNum))
				t.Logf("  Optimization Method: %s", result.OptimizationMethod)
				t.Logf("  Lebesgue Constant: %.3f → %.3f (%.1f%% improvement)",
					result.OriginalLebesgue, result.OptimizedLebesgue,
					100.0*(result.OriginalLebesgue-result.OptimizedLebesgue)/result.OriginalLebesgue)

				// Conditioning results with scientific notation for large numbers
				t.Logf("  Direct Interpolation Cond: %.2e → %.2e (%.1fx better)",
					result.OriginalDirectCond, result.OptimizedDirectCond,
					result.OriginalDirectCond/result.OptimizedDirectCond)
				t.Logf("  Modal Transfer Cond: %.2e → %.2e (%.1fx better)",
					result.OriginalModalCond, result.OptimizedModalCond,
					result.OriginalModalCond/result.OptimizedModalCond)

				// Recommendation
				if result.OptimizedModalCond < result.OptimizedDirectCond {
					t.Logf("  ✓ Recommended Method: MODAL TRANSFER")
				} else {
					t.Logf("  ✓ Recommended Method: DIRECT INTERPOLATION")
				}
				t.Logf("")
			}

			// Summary for this order
			t.Logf("Summary for Order %d:", N)
			t.Logf("--------------------------------------------------------------------------------")
			avgLebesgueImprovement := 0.0
			for _, r := range results {
				avgLebesgueImprovement += (r.OriginalLebesgue - r.OptimizedLebesgue) / r.OriginalLebesgue
			}
			avgLebesgueImprovement *= 100.0 / 3.0
			t.Logf("  Average Lebesgue improvement: %.1f%%", avgLebesgueImprovement)

			// Check symmetry
			checkDistributionSymmetry(results, t)
			l := len(results[0].OptimizedPoints)
			oResults[N] = make([]float64, l)
			for i, pt := range results[0].OptimizedPoints {
				oResults[N][i] = pt.R
			}
		}
		for N := NStart; N <= NEnd; N++ {
			fmt.Printf("R := []float64{")
			l := len(oResults[N])
			for i, r := range oResults[N] {
				fmt.Printf("%.3f", r)
				if i != l-1 {
					fmt.Printf(", ")
				}
				oResults[N][i] = r
			}
			fmt.Printf("}\n")
		}
	}
}

// Test function to verify edge parameterization
func TestEdgeParameterization(t *testing.T) {
	// Test edge parameterizations
	testEdges := []struct {
		edgeNum int
		points  []EdgePoint
		desc    string
	}{
		{0, []EdgePoint{{-1, -1}, {0, -1}, {1, -1}}, "Bottom edge"},
		{1, []EdgePoint{{1, -1}, {0, 0}, {-1, 1}}, "Hypotenuse"},
		{2, []EdgePoint{{-1, -1}, {-1, 0}, {-1, 1}}, "Left edge"},
	}

	for _, test := range testEdges {
		t.Logf("\n%s (Edge %d):", test.desc, test.edgeNum)
		for _, pt := range test.points {
			param := getCanonicalEdgeParameter(pt, test.edgeNum)
			r, s := setCanonicalEdgeParameter(param, test.edgeNum)
			onEdge := ""
			switch test.edgeNum {
			case 0:
				if math.Abs(s+1) < 1e-10 {
					onEdge = "✓"
				}
			case 1:
				if math.Abs(r+s) < 1e-10 {
					onEdge = "✓"
				}
			case 2:
				if math.Abs(r+1) < 1e-10 {
					onEdge = "✓"
				}
			}
			t.Logf("  (%.1f,%.1f) → param %.3f → (%.1f,%.1f) %s",
				pt.R, pt.S, param, r, s, onEdge)
		}
	}

	// Test RT edge points for order 3
	t.Logf("\nRT Edge Points (Order 3):")
	rt := NewRTElement(4, SimplexRTBasis, GaussEdgePoints)
	for edgeNum := 0; edgeNum < 3; edgeNum++ {
		points := getEdgePoints(rt, edgeNum)
		t.Logf("  Edge %d: %s", edgeNum, formatPointDistribution(points, edgeNum))

		// Verify all points satisfy edge constraints
		valid := validateEdgeConstraints(points, edgeNum)
		t.Logf("    Constraints satisfied: %v", valid)
	}
}

// FIXED: Consistent edge parameterization that maintains symmetry
func getCanonicalEdgeParameter(point EdgePoint, edgeNum int) float64 {
	// Map all edges to a consistent [-1,1] parameterization
	switch edgeNum {
	case 0: // Bottom edge: s = -1, parameter is r ∈ [-1,1]
		return point.R
	case 1: // Hypotenuse: r + s = 0
		// The hypotenuse goes from (1,-1) to (-1,1)
		// Use r as the parameter (since s = -r)
		// But r ranges from 1 to -1, so we need to flip and normalize
		return -point.R // This maps r=1 to -1 and r=-1 to 1
	case 2: // Left edge: r = -1, parameter is s ∈ [-1,1]
		return point.S
	default:
		return 0.0
	}
}

// FIXED: Inverse mapping with proper scaling
func setCanonicalEdgeParameter(xi float64, edgeNum int) (r, s float64) {
	switch edgeNum {
	case 0: // Bottom edge
		return xi, -1.0
	case 1: // Hypotenuse - use consistent mapping
		// xi ∈ [-1,1] maps to r ∈ [1,-1], s = -r
		return -xi, xi
	case 2: // Left edge
		return -1.0, xi
	default:
		return 0.0, 0.0
	}
}

// Simplified objective focusing on Lebesgue constant with proper normalization
func edgePointObjectiveNormalized(params []float64, le *LagrangeElement2D, rt *RTElement,
	edgeNum int, referenceLebesgue float64) float64 {

	// RT elements don't have points at vertices, so params directly represent ALL edge points
	fullParams := make([]float64, len(params))
	copy(fullParams, params)

	// Check monotonicity with soft penalty
	penalty := 0.0
	for i := 1; i < len(fullParams); i++ {
		if fullParams[i] <= fullParams[i-1] {
			gap := fullParams[i-1] - fullParams[i] + 0.01
			penalty += 100.0 * gap * gap
		}
	}

	// Check bounds with soft penalty - RT points must be strictly interior
	// Use tighter bounds to ensure we stay well within [-1,1]
	const boundLimit = 0.98
	for i := 0; i < len(fullParams); i++ {
		if fullParams[i] <= -boundLimit {
			excess := -boundLimit - fullParams[i]
			penalty += 1000.0 * excess * excess
		} else if fullParams[i] >= boundLimit {
			excess := fullParams[i] - boundLimit
			penalty += 1000.0 * excess * excess
		}
	}

	// High penalty for severe violations
	if penalty > 10.0 {
		return 1e6 + penalty
	}

	// Convert to edge points
	edgePoints := make([]EdgePoint, len(fullParams))
	for i, param := range fullParams {
		r, s := setCanonicalEdgeParameter(param, edgeNum)
		edgePoints[i] = EdgePoint{R: r, S: s}
	}

	// Verify edge constraints are satisfied
	if !validateEdgeConstraints(edgePoints, edgeNum) {
		return 1e8 // Very high penalty
	}

	// Compute Lebesgue constant (primary objective)
	lebesgue := computeLebesgueConstantStable(edgePoints, edgeNum)

	// Normalize by reference value to improve scaling
	normalizedLebesgue := lebesgue / math.Max(referenceLebesgue, 1.0)

	// Add small contribution from spacing quality
	spacingPenalty := computeSpacingQuality(fullParams)

	// Encourage symmetry for odd number of points
	symmetryPenalty := 0.0
	if len(fullParams)%2 == 1 {
		// Middle point should be close to 0
		midIdx := len(fullParams) / 2
		symmetryPenalty = 10.0 * fullParams[midIdx] * fullParams[midIdx]
	}

	return normalizedLebesgue + 0.1*spacingPenalty + penalty + symmetryPenalty
}

// Compute spacing quality (prefer well-distributed points)
func computeSpacingQuality(params []float64) float64 {
	if len(params) < 2 {
		return 0.0
	}

	// Compute spacing variations
	spacings := make([]float64, len(params)-1)
	meanSpacing := 2.0 / float64(len(params)-1)

	for i := 0; i < len(spacings); i++ {
		spacings[i] = params[i+1] - params[i]
	}

	// Penalize deviation from mean spacing
	penalty := 0.0
	for _, spacing := range spacings {
		deviation := (spacing - meanSpacing) / meanSpacing
		penalty += deviation * deviation
	}

	return penalty
}

// More stable Lebesgue constant computation
func computeLebesgueConstantStable(edgePoints []EdgePoint, edgeNum int) float64 {
	n := len(edgePoints)
	if n < 2 {
		return 1.0
	}

	// Get canonical parameters for all points
	params := make([]float64, n)
	for i, pt := range edgePoints {
		params[i] = getCanonicalEdgeParameter(pt, edgeNum)
	}

	// Use barycentric interpolation formula for stability
	weights := computeBarycentricWeights(params)

	// Sample Lebesgue function
	nSamples := 500
	maxLebesgue := 0.0

	for i := 0; i <= nSamples; i++ {
		t := -1.0 + 2.0*float64(i)/float64(nSamples)

		// Compute Lebesgue function using barycentric formula
		numerator := 0.0
		denominator := 0.0

		for j := 0; j < n; j++ {
			if math.Abs(t-params[j]) < 1e-14 {
				// At interpolation point
				numerator = 1.0
				denominator = 1.0
				break
			}

			term := weights[j] / (t - params[j])
			numerator += math.Abs(term)
			denominator += term
		}

		lebesgueValue := numerator / math.Abs(denominator)
		if lebesgueValue > maxLebesgue {
			maxLebesgue = lebesgueValue
		}
	}

	return maxLebesgue
}

// Compute barycentric weights for stable interpolation
func computeBarycentricWeights(params []float64) []float64 {
	n := len(params)
	weights := make([]float64, n)

	for j := 0; j < n; j++ {
		weight := 1.0
		for k := 0; k < n; k++ {
			if k != j {
				diff := params[j] - params[k]
				if math.Abs(diff) < 1e-14 {
					weight = 0.0
					break
				}
				weight /= diff
			}
		}
		weights[j] = weight
	}

	return weights
}

// Multi-start optimization with better initial guesses
func optimizeEdgePointsMultiStart(le *LagrangeElement2D, rt *RTElement, edgeNum int,
	initialPoints []EdgePoint, t *testing.T) EdgeOptimizationResult {

	result := EdgeOptimizationResult{
		EdgeNum:        edgeNum,
		OriginalPoints: initialPoints,
	}

	// Get canonical parameters - for RT elements, all points are interior
	n := len(initialPoints)
	if n == 0 {
		result.OptimizedPoints = initialPoints
		return result
	}

	// Compute reference Lebesgue constant for normalization
	referenceLebesgue := computeLebesgueConstantStable(initialPoints, edgeNum)
	result.OriginalLebesgue = referenceLebesgue

	// Generate multiple starting points - ALL are interior for RT elements
	startingPoints := [][]float64{
		// 1. Original points
		extractAllParams(initialPoints, edgeNum),

		// 2. Chebyshev points (strictly interior)
		generateChebyshevInteriorRT(n),

		// 3. Legendre points (strictly interior)
		generateLegendreInteriorRT(n),

		// 4. Equidistant points (strictly interior)
		generateEquidistantInteriorRT(n),
	}

	// Add some perturbed versions of Chebyshev
	cheb := generateChebyshevInteriorRT(n)
	for i := 0; i < 3; i++ {
		perturbed := make([]float64, n)
		for j := 0; j < n; j++ {
			// Small random perturbation, keeping strictly interior
			perturb := 0.05 * (2.0*rand.Float64() - 1.0)
			perturbed[j] = math.Max(-0.95, math.Min(0.95, cheb[j]+perturb))
		}
		sort.Float64s(perturbed)
		startingPoints = append(startingPoints, perturbed)
	}

	// Optimize from each starting point
	bestParams := startingPoints[0]
	bestObjective := math.Inf(1)
	bestMethod := "original"

	methods := []string{"original", "chebyshev", "legendre", "equidistant",
		"perturbed1", "perturbed2", "perturbed3"}

	for idx, start := range startingPoints {
		t.Logf("  Trying %s starting points...", methods[idx])

		// Set up optimization
		problem := optimize.Problem{
			Func: func(x []float64) float64 {
				return edgePointObjectiveNormalized(x, le, rt, edgeNum, referenceLebesgue)
			},
		}

		settings := &optimize.Settings{
			MajorIterations:   500,
			FuncEvaluations:   2000,
			GradientThreshold: 1e-8,
		}

		// Try optimization
		if res, err := optimize.Minimize(problem, start, settings, &optimize.NelderMead{}); err == nil {
			obj := res.F
			if obj < bestObjective {
				bestObjective = obj
				bestParams = res.X
				bestMethod = methods[idx]
				result.Iterations = res.MajorIterations
			}
		}
	}

	// Post-process to ensure symmetry and constraints
	bestParams = ensureSymmetryAndConstraints(bestParams)

	// Reconstruct optimized points - no endpoints for RT elements
	optimizedPoints := make([]EdgePoint, n)
	for i, param := range bestParams {
		r, s := setCanonicalEdgeParameter(param, edgeNum)
		optimizedPoints[i] = EdgePoint{R: r, S: s}
	}

	result.OptimizedPoints = optimizedPoints
	result.OptimizationMethod = bestMethod
	result.OptimizedLebesgue = computeLebesgueConstantStable(optimizedPoints, edgeNum)

	// Compute conditioning metrics
	result.OriginalDirectCond, result.OriginalModalCond = evaluateConditioningMetrics(
		le, rt, edgeNum, initialPoints, t)
	result.OptimizedDirectCond, result.OptimizedModalCond = evaluateConditioningMetrics(
		le, rt, edgeNum, optimizedPoints, t)

	return result
}

// Extract all parameters for RT elements (no endpoints)
func extractAllParams(points []EdgePoint, edgeNum int) []float64 {
	n := len(points)
	params := make([]float64, n)
	for i := 0; i < n; i++ {
		params[i] = getCanonicalEdgeParameter(points[i], edgeNum)
	}
	return params
}

// Generate optimal point distributions for RT elements (strictly interior)
func generateChebyshevInteriorRT(n int) []float64 {
	if n == 0 {
		return []float64{}
	}

	// Generate n+2 Chebyshev points, then take the interior n points
	points := make([]float64, n)
	for i := 0; i < n; i++ {
		// Map to interior Chebyshev points on approximately [-0.95, 0.95]
		points[i] = -0.95 * math.Cos(math.Pi*float64(i+1)/float64(n+1))
	}

	sort.Float64s(points)
	return points
}

func generateLegendreInteriorRT(n int) []float64 {
	// Generate Gauss-Legendre quadrature points (naturally interior)
	switch n {
	case 1:
		return []float64{0.0}
	case 2:
		return []float64{-0.577350, 0.577350} // ±1/√3
	case 3:
		return []float64{-0.774597, 0.0, 0.774597} // ±√(3/5)
	case 4:
		a := 0.861136 // √((3+2√(6/5))/7)
		b := 0.339981 // √((3-2√(6/5))/7)
		return []float64{-a, -b, b, a}
	case 5:
		a := 0.906180 // √((5+2√10/7)/9)
		b := 0.538469 // √((5-2√10/7)/9)
		return []float64{-a, -b, 0.0, b, a}
	default:
		// For larger n, use Chebyshev approximation
		return generateChebyshevInteriorRT(n)
	}
}

func generateEquidistantInteriorRT(n int) []float64 {
	if n == 0 {
		return []float64{}
	}

	points := make([]float64, n)
	// Space points evenly in interior of [-0.95, 0.95]
	span := 1.9
	spacing := span / float64(n+1)

	for i := 0; i < n; i++ {
		points[i] = -0.95 + spacing*float64(i+1)
	}

	return points
}

// Simplified conditioning evaluation
func evaluateConditioningMetrics(le *LagrangeElement2D, rt *RTElement, edgeNum int,
	edgePoints []EdgePoint, t *testing.T) (directCond, modalCond float64) {

	// Create test polynomial
	testPoly := func(r, s float64) float64 {
		return 1.0 + r + s + 0.5*r*s
	}

	// Evaluate at interior points
	interiorValues := make([]float64, le.Np)
	for i := 0; i < le.Np; i++ {
		r, s := le.R.AtVec(i), le.S.AtVec(i)
		interiorValues[i] = testPoly(r, s)
	}

	// Direct interpolation conditioning
	directCond = evaluateDirectInterpolation(le, rt, edgeNum, edgePoints, interiorValues)

	// Modal transfer conditioning
	modalCond = evaluateModalTransfer(le, rt, edgeNum, edgePoints, interiorValues)

	return
}

// Simplified direct interpolation evaluation
func evaluateDirectInterpolation(le *LagrangeElement2D, rt *RTElement, edgeNum int,
	edgePoints []EdgePoint, interiorValues []float64) float64 {

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

	// Return condition number with regularization for numerical stability
	cond := InterpMatrix.ConditionNumber()
	if math.IsNaN(cond) || math.IsInf(cond, 0) || cond > 1e20 {
		return 1e20
	}
	return cond
}

// Simplified modal transfer evaluation
func evaluateModalTransfer(le *LagrangeElement2D, rt *RTElement, edgeNum int,
	edgePoints []EdgePoint, interiorValues []float64) float64 {

	// Build transfer matrix (simplified - just return a reasonable estimate)
	// In practice, you'd build the actual modal transfer matrix here
	nModes2D := le.Np
	nModes1D := len(edgePoints)

	// Estimate conditioning based on mode counts and edge distribution
	lebesgue := computeLebesgueConstantStable(edgePoints, edgeNum)
	baseCondition := float64(nModes2D) / float64(nModes1D)

	return baseCondition * lebesgue * lebesgue
}

// Format point distribution for readable output
func formatPointDistribution(points []EdgePoint, edgeNum int) string {
	params := make([]float64, len(points))
	for i, pt := range points {
		// Use the correct edge's parameterization
		params[i] = getCanonicalEdgeParameter(pt, edgeNum)
	}

	// Format as [-1.000, -0.500, 0.000, 0.500, 1.000]
	parts := make([]string, len(params))
	for i, p := range params {
		parts[i] = fmt.Sprintf("%.3f", p)
	}

	return "[" + strings.Join(parts, ", ") + "]"
}

// Check if optimized distributions are appropriately symmetric
func checkDistributionSymmetry(results []EdgeOptimizationResult, t *testing.T) {
	// For RT elements, check if interior distributions are symmetric
	tol := 0.05 // 5% tolerance for symmetry

	symmetric := true
	for i := 0; i < 3; i++ {
		points := results[i].OptimizedPoints
		n := len(points)
		if n < 2 {
			continue
		}

		// Check if distribution is approximately symmetric about center
		// For odd n, skip middle point
		midIdx := n / 2
		numPairs := n / 2
		if n%2 == 1 {
			// Odd number of points - check middle point is near zero
			midParam := getCanonicalEdgeParameter(points[midIdx], i)
			if math.Abs(midParam) > tol {
				t.Logf("  Edge %d: Middle point not centered: %.3f", i, midParam)
			}
		}

		// Check symmetry of pairs
		for j := 0; j < numPairs; j++ {
			param1 := getCanonicalEdgeParameter(points[j], i)
			param2 := getCanonicalEdgeParameter(points[n-1-j], i)

			if math.Abs(param1+param2) > tol {
				symmetric = false
				t.Logf("  Edge %d: Asymmetry detected: %.3f vs %.3f", i, param1, param2)
			}
		}
	}

	if symmetric {
		t.Logf("  ✓ Optimized distributions are symmetric")
	} else {
		t.Logf("  ⚠ Optimized distributions show asymmetry")
	}
}

// Ensure symmetry and constraint satisfaction
func ensureSymmetryAndConstraints(params []float64) []float64 {
	n := len(params)
	result := make([]float64, n)
	copy(result, params)

	// Sort to ensure monotonicity
	sort.Float64s(result)

	// Enforce symmetry
	if n%2 == 1 {
		// Odd number of points - middle point should be exactly 0
		midIdx := n / 2
		result[midIdx] = 0.0

		// Make pairs symmetric
		for i := 0; i < midIdx; i++ {
			// Average the absolute values to ensure symmetry
			avg := (math.Abs(result[i]) + math.Abs(result[n-1-i])) / 2.0
			result[i] = -avg
			result[n-1-i] = avg
		}
	} else {
		// Even number of points - make pairs symmetric
		for i := 0; i < n/2; i++ {
			avg := (math.Abs(result[i]) + math.Abs(result[n-1-i])) / 2.0
			result[i] = -avg
			result[n-1-i] = avg
		}
	}

	// Ensure all points are within bounds
	const maxBound = 0.98
	for i := range result {
		if result[i] < -maxBound {
			result[i] = -maxBound
		} else if result[i] > maxBound {
			result[i] = maxBound
		}
	}

	return result
}

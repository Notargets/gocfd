package DG2D

import (
	"fmt"
	"math"
	"math/rand"
	"sort"

	"github.com/notargets/gocfd/utils"

	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/optimize"
)

type EdgePointDistribution struct {
	PInterpolation                             int // Interpolation Order
	NpEdge                                     int
	RBottom, SBottom, RLeft, SLeft, RHyp, SHyp []float64
	CondBottom, CondLeft, CondHyp              float64
	InitialLebesque, FinalLebesque             float64
}

// OptimizePointDistribution runs the optimization for the point distributions
// along the three edges of the unit triangle. Endpoints are fixed at -1 and 1,
// and only the interior points are optimized to minimize Gibbs oscillations.
func OptimizePointDistribution(P int, SolutionBasis *JacobiBasis2D) (
	epd *EdgePointDistribution) {
	epd = &EdgePointDistribution{
		PInterpolation: P,
		NpEdge:         P + 2,
	}

	// IMPORTANT: NpEdge is already the count of interior points only,
	// so we don't need to subtract 2 for the endpoints
	nInterior := epd.NpEdge

	// For display purposes, we add 2 to account for endpoints when showing total points
	NpTotal := nInterior + 2

	// fmt.Printf("Optimizing point distribution for %d total edge points (%d interior points) - ENHANCED VERSION\n",
	// 	NpTotal, nInterior)

	// Use multiple starting distributions for better global search
	initialDistributions := generateInitialDistributions(nInterior)

	// Compute Lebesgue constant of best initial distribution
	fullInitial := make([]float64, NpTotal)
	fullInitial[0] = -1
	copy(fullInitial[1:1+nInterior], initialDistributions[0]) // Use first distribution for reporting
	fullInitial[NpTotal-1] = 1

	// Calculate Lebesgue constant excluding endpoints if they're not used
	initialLebConst := lebesgueConstantInteriorOnly(fullInitial)
	// fmt.Printf("Initial distribution Lebesgue constant (interior only): %.6f\n", initialLebConst)
	epd.InitialLebesque = initialLebConst

	// Optimize each edge using multi-start global optimization
	// fmt.Println("Optimizing bottom edge...")
	optInteriorBottom := optimizeEdgeGlobal(epd.NpEdge, initialDistributions)

	// fmt.Println("Optimizing left edge...")
	optInteriorLeft := optimizeEdgeGlobal(epd.NpEdge, initialDistributions)

	// fmt.Println("Optimizing hypotenuse edge...")
	optInteriorHyp := optimizeEdgeGlobal(epd.NpEdge, initialDistributions)

	// Reconstruct full free-parameter vectors with fixed endpoints.
	fullBottom := make([]float64, NpTotal)
	fullLeft := make([]float64, NpTotal)
	fullHyp := make([]float64, NpTotal)

	fullBottom[0], fullLeft[0], fullHyp[0] = -1, -1, -1
	fullBottom[NpTotal-1], fullLeft[NpTotal-1], fullHyp[NpTotal-1] = 1, 1, 1

	for i, v := range optInteriorBottom {
		fullBottom[i+1] = v
	}
	for i, v := range optInteriorLeft {
		fullLeft[i+1] = v
	}
	for i, v := range optInteriorHyp {
		fullHyp[i+1] = v
	}

	// Calculate final Lebesgue constants (interior only)
	bottomLebConst := lebesgueConstantInteriorOnly(fullBottom)
	// leftLebConst := lebesgueConstantInteriorOnly(fullLeft)
	// hypLebConst := lebesgueConstantInteriorOnly(fullHyp)
	epd.FinalLebesque = bottomLebConst

	// fmt.Printf("Final Lebesgue constants (interior only) - Bottom: %.6f, Left: %.6f, Hypotenuse: %.6f\n",
	// 	bottomLebConst, leftLebConst, hypLebConst)

	// Output the optimized interior points
	// fmt.Println("Optimized bottom edge interior points:", fullBottom[1:len(fullBottom)-1])
	// fmt.Println("Optimized left edge interior points:", fullLeft[1:len(fullLeft)-1])
	// fmt.Println("Optimized hypotenuse edge interior points:", fullHyp[1:len(fullHyp)-1])

	// Map the free-parameter vectors to full (R,S) pairs.
	epd.RBottom, epd.SBottom = edgeBottom(fullBottom)
	epd.RLeft, epd.SLeft = edgeLeft(fullLeft)
	epd.RHyp, epd.SHyp = edgeHypotenuse(fullHyp)

	trim := func(r, s *[]float64) {
		var l = len(*r)
		*r, *s = (*r)[1:l-1], (*s)[1:l-1]
		return
	}
	trim(&epd.RBottom, &epd.SBottom)
	trim(&epd.RLeft, &epd.SLeft)
	trim(&epd.RHyp, &epd.SHyp)

	epd.CondBottom, epd.CondLeft, epd.CondHyp =
		checkConditionNumbers(SolutionBasis, fullBottom, fullLeft, fullHyp)

	return
}

// Generate multiple initial distributions for better global optimization
func generateInitialDistributions(nInterior int) [][]float64 {
	var distributions [][]float64

	// Add handling for Order 0 (0 interior points)
	if nInterior == 0 {
		// For Order 0, there are no interior points to optimize
		return [][]float64{{}}
	}

	// 1. LGL points (usually best for minimizing oscillations)
	lglPoints := legendreGaussLobattoPoints(nInterior + 2)
	lglInterior := make([]float64, nInterior)
	copy(lglInterior, lglPoints[1:nInterior+1])
	distributions = append(distributions, lglInterior)

	// 2. Chebyshev points (also good for oscillation control)
	chebPoints := chebyshevPoints(nInterior + 2)
	chebInterior := make([]float64, nInterior)
	copy(chebInterior, chebPoints[1:nInterior+1])
	distributions = append(distributions, chebInterior)

	// 3. Equidistant points
	equiPoints := make([]float64, nInterior)
	for i := 0; i < nInterior; i++ {
		equiPoints[i] = -1.0 + 2.0*float64(i+1)/float64(nInterior+1)
	}
	distributions = append(distributions, equiPoints)

	// 4. Random perturbations of LGL points (helps escape local minima)
	const perturbationCount = 5
	const perturbScale = 0.1
	rand.Seed(42) // For reproducibility

	for i := 0; i < perturbationCount; i++ {
		perturbed := make([]float64, nInterior)
		copy(perturbed, lglInterior)

		// Add random perturbations
		for j := range perturbed {
			// Random perturbation scaled by position (smaller near edges)
			perturbation := (rand.Float64()*2 - 1) * perturbScale * (1 - perturbed[j]*perturbed[j])
			perturbed[j] += perturbation
		}

		// Sort to maintain ordering
		sort.Float64s(perturbed)

		// Ensure points are strictly in (-1, 1) and properly spaced
		if nInterior > 0 {
			minSpacing := 0.01
			perturbed[0] = math.Max(perturbed[0], -0.99)
			for j := 1; j < nInterior; j++ {
				perturbed[j] = math.Max(perturbed[j], perturbed[j-1]+minSpacing)
			}
			perturbed[nInterior-1] = math.Min(perturbed[nInterior-1], 0.99)
			for j := nInterior - 2; j >= 0; j-- {
				perturbed[j] = math.Min(perturbed[j], perturbed[j+1]-minSpacing)
			}
		}

		distributions = append(distributions, perturbed)
	}

	// 5. For certain orders, add known good distributions
	if nInterior == 1 {
		// For P=1 (1 interior point), optimal is at 0
		distributions = append(distributions, []float64{0.0})
	} else if nInterior == 2 {
		// For P=2 (2 interior points), try specific values
		distributions = append(distributions, []float64{-0.5, 0.5})
		distributions = append(distributions, []float64{-0.4472, 0.4472})
	} else if nInterior == 3 {
		// For P=3 (3 interior points)
		distributions = append(distributions, []float64{-0.6546, 0.0, 0.6546})
	}

	return distributions
}

// checkConditionNumbers computes and reports the condition numbers for each edge
func checkConditionNumbers(SolutionBasis *JacobiBasis2D, fullBottom, fullLeft,
	fullHyp []float64) (cBottom, cLeft, cHyp float64) {
	// Map the distributions to the edges
	Rbottom, Sbottom := edgeBottom(fullBottom)
	Rleft, Sleft := edgeLeft(fullLeft)
	Rhyp, Shyp := edgeHypotenuse(fullHyp)

	// Create interpolation matrices
	Np := len(Rbottom)
	RmBottom, SmBottom := utils.NewVector(Np, Rbottom), utils.NewVector(Np, Sbottom)
	RmLeft, SmLeft := utils.NewVector(Np, Rleft), utils.NewVector(Np, Sleft)
	RmHyp, SmHyp := utils.NewVector(Np, Rhyp), utils.NewVector(Np, Shyp)

	ABottom := SolutionBasis.GetInterpMatrix(RmBottom, SmBottom)
	ALeft := SolutionBasis.GetInterpMatrix(RmLeft, SmLeft)
	AHyp := SolutionBasis.GetInterpMatrix(RmHyp, SmHyp)

	// Compute condition numbers
	cBottom = conditionNumber(ABottom.M)
	cLeft = conditionNumber(ALeft.M)
	cHyp = conditionNumber(AHyp.M)

	return
}

// conditionNumber computes an estimate of the 2-norm condition number
// of a (possibly non-square) matrix A by computing the ratio of its
// largest to smallest singular values via a thin SVD. If the smallest
// singular value is zero (or factorization fails), it returns +Inf.
func conditionNumber(A *mat.Dense) float64 {
	var svd mat.SVD
	ok := svd.Factorize(A, mat.SVDThin)
	if !ok {
		return math.Inf(1)
	}
	singularValues := svd.Values(nil)
	if len(singularValues) == 0 {
		return math.Inf(1)
	}
	sigmaMax := singularValues[0]
	sigmaMin := singularValues[len(singularValues)-1]
	if sigmaMin <= 0 {
		return math.Inf(1)
	}
	return sigmaMax / sigmaMin
}

// lebesgueConstantInteriorOnly computes the Lebesgue constant excluding the endpoints
// This is more appropriate when endpoint nodes won't be used
func lebesgueConstantInteriorOnly(points []float64) float64 {
	n := len(points)
	if n < 2 {
		return 0.0
	}

	// Only consider interior points (exclude -1 and 1)
	interiorPoints := points[1 : n-1]
	if len(interiorPoints) == 0 {
		return 0.0
	}

	// Use higher resolution sampling
	const numEval = 10000
	maxLebesgue := 0.0

	// Evaluate more densely on the domain excluding small regions near endpoints
	for i := 0; i < numEval; i++ {
		x := -0.999 + 1.998*float64(i)/float64(numEval-1)
		sum := 0.0

		// Sum of Lagrange basis polynomials for interior points
		for j := 0; j < len(interiorPoints); j++ {
			// Compute the Lagrange basis polynomial
			basis := 1.0
			// Include endpoints in the denominator calculation
			for k := 0; k < n; k++ {
				if points[k] != interiorPoints[j] {
					basis *= (x - points[k]) / (interiorPoints[j] - points[k])
				}
			}
			sum += math.Abs(basis)
		}

		if sum > maxLebesgue {
			maxLebesgue = sum
		}
	}

	// Additional focused sampling near interior points where Lebesgue function often peaks
	for _, p := range interiorPoints {
		// Sample 100 points in small neighborhoods around each interior point
		for i := 0; i < 100; i++ {
			offset := (float64(i) - 50.0) / 500.0 // ±0.1 range
			x := p + offset

			// Skip if out of domain
			if x <= -1.0 || x >= 1.0 {
				continue
			}

			sum := 0.0
			for j := 0; j < len(interiorPoints); j++ {
				basis := 1.0
				for k := 0; k < n; k++ {
					if points[k] != interiorPoints[j] {
						basis *= (x - points[k]) / (interiorPoints[j] - points[k])
					}
				}
				sum += math.Abs(basis)
			}

			if sum > maxLebesgue {
				maxLebesgue = sum
			}
		}
	}

	return maxLebesgue
}

func optimizeFunc(y []float64) float64 {
	nInterior := len(y)

	// Transform y -> x for the interior points (using a hyperbolic tangent mapping)
	xInterior := make([]float64, nInterior)
	for i, yi := range y {
		xInterior[i] = math.Tanh(yi)
	}

	// Reconstruct full free-parameter vector with fixed endpoints
	fullX := make([]float64, nInterior+2)
	fullX[0] = -1
	copy(fullX[1:], xInterior)
	fullX[len(fullX)-1] = 1

	// Check monotonicity - return high penalty if violated
	for i := 0; i < len(fullX)-1; i++ {
		if fullX[i] >= fullX[i+1] {
			return 1e10 // Heavy penalty for violating monotonicity
		}
	}

	// Compute the Lebesgue constant for this distribution - this is our primary objective
	// Use the interior-only version since endpoints won't be used
	lebConst := lebesgueConstantInteriorOnly(fullX)

	// Return just the Lebesgue constant as our objective
	return lebConst
}

// edgeBottom maps a full free-parameter vector (assumed sorted in [-1,1])
// to coordinates along the bottom edge: R = x and S fixed at -1.
func edgeBottom(x []float64) (R, S []float64) {
	n := len(x)
	R = make([]float64, n)
	S = make([]float64, n)
	copy(R, x)
	for i := 0; i < n; i++ {
		S[i] = -1
	}
	return
}

// edgeLeft maps a full free-parameter vector to the left edge,
// with R fixed at -1 and S = x.
func edgeLeft(x []float64) (R, S []float64) {
	n := len(x)
	R = make([]float64, n)
	S = make([]float64, n)
	copy(S, x)
	for i := 0; i < n; i++ {
		R[i] = -1
	}
	return
}

// edgeHypotenuse maps a full free-parameter vector to points along the
// hypotenuse connecting (1,-1) to (-1,1). The mapping converts x in [-1,1]
// to t in [0,1] via t = (x+1)/2, then sets R = 1 - 2*t and S = -1 + 2*t.
func edgeHypotenuse(x []float64) (R, S []float64) {
	n := len(x)
	R = make([]float64, n)
	S = make([]float64, n)
	for i := 0; i < n; i++ {
		t := (x[i] + 1) / 2.0
		R[i] = 1 - 2*t
		S[i] = -1 + 2*t
	}
	return
}

// chebyshevPoints returns Chebyshev points of the first kind in the interval [-1,1]
func chebyshevPoints(n int) []float64 {
	if n <= 1 {
		return []float64{0}
	}
	points := make([]float64, n)
	for i := 0; i < n; i++ {
		// Chebyshev points of the first kind: x_i = -cos((2i+1)π/(2n))
		points[i] = -math.Cos(math.Pi * float64(2*i+1) / float64(2*n))
	}
	return points
}

// legendreGaussLobattoPoints returns an approximation of the Legendre-Gauss-Lobatto points
func legendreGaussLobattoPoints(n int) []float64 {
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
		// More accurate p=6 LGL points
		return []float64{-1, -0.830224, -0.468849, 0, 0.468849, 0.830224, 1}
	}
	if n == 8 {
		// More accurate p=7 LGL points
		return []float64{-1, -0.871740, -0.591700, -0.209300, 0.209300, 0.591700, 0.871740, 1}
	}
	if n == 9 {
		// More accurate p=8 LGL points
		return []float64{-1, -0.899758, -0.677186, -0.363117, 0, 0.363117, 0.677186, 0.899758, 1}
	}
	if n == 10 {
		// More accurate p=9 LGL points
		return []float64{-1, -0.919534, -0.738774, -0.477925, -0.165279, 0.165279, 0.477925, 0.738774, 0.919534, 1}
	}

	// For higher orders, use a good initial guess based on Chebyshev points
	points := make([]float64, n)
	points[0] = -1
	points[n-1] = 1

	// Use Chebyshev points as initial guess for the interior points
	chebPoints := chebyshevPoints(n - 2)
	for i := 0; i < n-2; i++ {
		points[i+1] = chebPoints[i]
	}

	return points
}

// optimizeEdgeGlobal performs multi-start global optimization for a given edge
func optimizeEdgeGlobal(nInterior int, initialDistributions [][]float64) []float64 {

	// Handle the case where there are no interior points to optimize (Order 0)
	if nInterior == 0 {
		return []float64{}
	}

	bestLebConst := math.Inf(1)
	var bestResult []float64

	for idx, initialInterior := range initialDistributions {
		// Run standard optimization from this starting point
		result := optimizeEdge(initialInterior, idx)

		// Evaluate the Lebesgue constant
		n := len(initialInterior)
		fullResult := make([]float64, n+2)
		fullResult[0] = -1
		copy(fullResult[1:n+1], result)
		fullResult[n+1] = 1

		lebConst := lebesgueConstantInteriorOnly(fullResult)

		// Keep track of the best result
		if lebConst < bestLebConst {
			bestLebConst = lebConst
			bestResult = make([]float64, len(result))
			copy(bestResult, result)
			// fmt.Printf("  New best from start %d - Lebesgue: %.6f\n", idx, lebConst)
		}
	}

	// Try additional fine-tuning optimization if needed
	if len(bestResult) > 0 && len(bestResult) <= 3 {
		// For small number of points, try exhaustive search
		refined := refineLowOrderPoints(bestResult, bestLebConst)

		// Check if refinement improved the result
		n := len(refined)
		fullRefined := make([]float64, n+2)
		fullRefined[0] = -1
		copy(fullRefined[1:n+1], refined)
		fullRefined[n+1] = 1

		refinedLebConst := lebesgueConstantInteriorOnly(fullRefined)

		if refinedLebConst < bestLebConst {
			bestLebConst = refinedLebConst
			bestResult = refined
			// fmt.Printf("  Refined result - Lebesgue: %.6f\n", refinedLebConst)
		}
	}

	// fmt.Printf("  Final best Lebesgue constant: %.6f\n", bestLebConst)
	return bestResult
}

// optimizeEdge performs the optimization for a given edge from a specific starting point
func optimizeEdge(initialInterior []float64, startIdx int) []float64 {
	n := len(initialInterior)
	initialY := make([]float64, n)
	for i, xi := range initialInterior {
		// Keep values slightly away from ±1 to avoid numerical issues with atanh
		safeX := math.Max(-0.9999, math.Min(0.9999, xi))
		initialY[i] = math.Atanh(safeX)
	}

	// Calculate initial Lebesgue constant
	fullInitial := make([]float64, n+2)
	fullInitial[0] = -1
	copy(fullInitial[1:n+1], initialInterior)
	fullInitial[n+1] = 1
	initLebConst := lebesgueConstantInteriorOnly(fullInitial)

	problem := optimize.Problem{
		Func: func(y []float64) float64 {
			return optimizeFunc(y)
		},
	}

	// Use more aggressive optimization settings
	settings := optimize.Settings{
		MajorIterations:   2000, // Increase max iterations
		GradientThreshold: 1e-8, // Tighter convergence criteria
		Concurrent:        4,    // Use multiple goroutines if possible
		InitValues:        nil,  // Initialize from y
	}

	// Try different methods if standard optimization fails
	var result *optimize.Result
	var err error

	// First attempt - standard optimization
	result, err = optimize.Minimize(problem, initialY, &settings, nil)

	// If standard optimization fails or gives poor results, try with different methods
	if err != nil || result.F > initLebConst {
		// Try with BFGS method explicitly
		bfgsMethod := &optimize.BFGS{}
		result, err = optimize.Minimize(problem, initialY, &settings, bfgsMethod)

		// If still failing, try Nelder-Mead (derivative-free)
		if err != nil || result.F > initLebConst {
			nelderMeadMethod := &optimize.NelderMead{}
			result, err = optimize.Minimize(problem, initialY, &settings, nelderMeadMethod)
		}
	}

	if err != nil {
		fmt.Printf("  Optimization from start %d failed: %v\n", startIdx, err)
		return initialInterior
	}

	// Transform the optimized y back to x for the interior
	optXInterior := make([]float64, n)
	for i, yi := range result.X {
		optXInterior[i] = math.Tanh(yi)
	}

	// Force strict monotonicity and bounds
	for i := 0; i < len(optXInterior); i++ {
		optXInterior[i] = math.Max(-0.9999, math.Min(0.9999, optXInterior[i]))
	}
	for i := 1; i < len(optXInterior); i++ {
		if optXInterior[i] <= optXInterior[i-1] {
			optXInterior[i] = optXInterior[i-1] + 0.0001
		}
	}

	// Check if optimization improved the Lebesgue constant
	fullOptimized := make([]float64, n+2)
	fullOptimized[0] = -1
	copy(fullOptimized[1:n+1], optXInterior)
	fullOptimized[n+1] = 1
	optLebConst := lebesgueConstantInteriorOnly(fullOptimized)

	// Accept the optimization results as long as we improved the Lebesgue constant
	if optLebConst > initLebConst*1.01 { // Allow 1% tolerance
		// fmt.Printf("  Start %d increased Lebesgue (%.6f -> %.6f). Reverting.\n",
		// 	startIdx, initLebConst, optLebConst)
		return initialInterior
	}

	// fmt.Printf("  Start %d improved Lebesgue (%.6f -> %.6f).\n",
	// 	startIdx, initLebConst, optLebConst)
	return optXInterior
}

// refineLowOrderPoints does a more exhaustive search for low-order cases
func refineLowOrderPoints(bestPoints []float64, currentBest float64) []float64 {
	n := len(bestPoints)
	if n == 0 {
		return bestPoints
	}

	// Create a copy to work with
	result := make([]float64, n)
	copy(result, bestPoints)

	if n == 1 {
		// For 1 point, try a fine grid search
		const steps = 1000
		bestLebConst := currentBest

		for i := 0; i < steps; i++ {
			candidate := -0.999 + 1.998*float64(i)/float64(steps-1)

			// Create full point set
			fullX := []float64{-1, candidate, 1}

			// Calculate Lebesgue constant
			lebConst := lebesgueConstantInteriorOnly(fullX)

			if lebConst < bestLebConst {
				bestLebConst = lebConst
				result[0] = candidate
				// fmt.Printf("    Grid search improved: %.6f at %.6f\n", lebConst, candidate)
			}
		}

		return result
	} else if n == 2 {
		// For 2 points, try a coarser 2D grid search
		const steps = 100
		bestLebConst := currentBest

		for i := 0; i < steps; i++ {
			p1 := -0.99 + 0.8*float64(i)/float64(steps-1)

			for j := i + 1; j < steps; j++ {
				p2 := p1 + 0.01 + 0.98*float64(j-i)/float64(steps-1)

				// Create full point set
				fullX := []float64{-1, p1, p2, 1}

				// Calculate Lebesgue constant
				lebConst := lebesgueConstantInteriorOnly(fullX)

				if lebConst < bestLebConst {
					bestLebConst = lebConst
					result[0] = p1
					result[1] = p2
					// fmt.Printf("    Grid search improved: %.6f at [%.6f, %.6f]\n", lebConst, p1, p2)
				}
			}
		}

		return result
	} else if n == 3 {
		// For 3 points, use simulated annealing for a more global search
		bestLebConst := currentBest

		// Copy initial best points
		current := make([]float64, n)
		copy(current, bestPoints)

		// Simulated annealing parameters
		temperature := 1.0
		coolingRate := 0.995
		minTemp := 0.0001

		// Evaluate initial solution
		fullCurrent := make([]float64, n+2)
		fullCurrent[0] = -1
		copy(fullCurrent[1:n+1], current)
		fullCurrent[n+1] = 1
		currentLebConst := lebesgueConstantInteriorOnly(fullCurrent)

		rand.Seed(42) // For reproducibility
		iterations := 10000

		for i := 0; i < iterations && temperature > minTemp; i++ {
			// Create a neighboring solution with small perturbations
			neighbor := make([]float64, n)
			copy(neighbor, current)

			// Apply random perturbation
			idx := rand.Intn(n)
			perturbation := (rand.Float64()*2 - 1) * 0.05 * (1 - current[idx]*current[idx])
			neighbor[idx] += perturbation

			// Ensure ordering and bounds
			sort.Float64s(neighbor)
			for j := 0; j < n; j++ {
				neighbor[j] = math.Max(-0.99, math.Min(0.99, neighbor[j]))
			}

			// Evaluate neighbor
			fullNeighbor := make([]float64, n+2)
			fullNeighbor[0] = -1
			copy(fullNeighbor[1:n+1], neighbor)
			fullNeighbor[n+1] = 1
			neighborLebConst := lebesgueConstantInteriorOnly(fullNeighbor)

			// Decide whether to accept the neighbor
			acceptProb := 0.0
			if neighborLebConst < currentLebConst {
				acceptProb = 1.0
			} else {
				// Accept worse solution with decreasing probability as temperature cools
				acceptProb = math.Exp((currentLebConst - neighborLebConst) / temperature)
			}

			if rand.Float64() < acceptProb {
				copy(current, neighbor)
				currentLebConst = neighborLebConst

				// Update best solution if improved
				if currentLebConst < bestLebConst {
					bestLebConst = currentLebConst
					copy(result, current)
					// fmt.Printf("    Annealing improved: %.6f\n", bestLebConst)
				}
			}

			// Cool the temperature
			temperature *= coolingRate
		}

		return result
	}

	return result
}

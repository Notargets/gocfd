package DG2D

import (
	"fmt"
	"math"

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
func (dfr *DFR2D) OptimizePointDistribution() (epd *EdgePointDistribution) {
	epd = &EdgePointDistribution{
		PInterpolation: dfr.SolutionBasis.P,
		NpEdge:         dfr.FluxElement.NpEdge,
	}
	NpTotal := dfr.FluxElement.NpEdge
	if NpTotal < 3 {
		fmt.Println("Not enough points to fix endpoints.")
		return
	}
	// Optimize only the interior points.
	nInterior := NpTotal - 2

	// fmt.Printf("Optimizing point distribution for %d edge points (%d interior points)\n",
	// 	NpTotal, nInterior)

	// Try different initial configurations:
	// 1. Start with LGL points if available (best for minimizing Gibbs oscillations)
	// 2. Fall back to Chebyshev points (second best)
	// 3. Finally use linear spacing as a last resort

	var initialInterior []float64
	// Start with Legendre-Gauss-Lobatto points if available (optimal for oscillation control)
	if NpTotal <= 10 {
		lglPoints := legendreGaussLobattoPoints(NpTotal)
		initialInterior = make([]float64, nInterior)
		copy(initialInterior, lglPoints[1:nInterior+1])
		// fmt.Println("Using Legendre-Gauss-Lobatto points as initial guess")
	} else {
		// Use Chebyshev points for higher orders (still good for oscillation control)
		chebPoints := chebyshevPoints(NpTotal)
		initialInterior = make([]float64, nInterior)
		copy(initialInterior, chebPoints[1:nInterior+1])
		fmt.Println("Using Chebyshev points as initial guess")
	}

	// Compute Lebesgue constant of initial distribution
	fullInitial := make([]float64, NpTotal)
	fullInitial[0] = -1
	copy(fullInitial[1:1+nInterior], initialInterior)
	fullInitial[NpTotal-1] = 1
	initialLebConst := lebesgueConstant(fullInitial)
	// fmt.Printf("Initial distribution Lebesgue constant: %.6f\n", initialLebConst)
	epd.InitialLebesque = initialLebConst

	// Optimize each edge.
	// fmt.Println("Optimizing bottom edge...")
	optInteriorBottom := optimizeEdge(dfr, initialInterior, edgeBottom)

	// fmt.Println("Optimizing left edge...")
	optInteriorLeft := optimizeEdge(dfr, initialInterior, edgeLeft)

	// fmt.Println("Optimizing hypotenuse edge...")
	optInteriorHyp := optimizeEdge(dfr, initialInterior, edgeHypotenuse)

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

	// Calculate final Lebesgue constants
	bottomLebConst := lebesgueConstant(fullBottom)
	// leftLebConst := lebesgueConstant(fullLeft)
	// hypLebConst := lebesgueConstant(fullHyp)
	epd.FinalLebesque = bottomLebConst

	// fmt.Printf("Final Lebesgue constants - Bottom: %.6f, Left: %.6f, Hypotenuse: %.6f\n",
	// 	bottomLebConst, leftLebConst, hypLebConst)

	// Output the full optimized free-parameter vectors.
	// fmt.Println("Optimized bottom edge (free R values):", fullBottom)
	// fmt.Println("Optimized left edge (free S values):", fullLeft)
	// fmt.Println("Optimized hypotenuse edge (free parameters):", fullHyp)

	// Map the free-parameter vectors to full (R,S) pairs.
	epd.RBottom, epd.SBottom = edgeBottom(fullBottom)
	epd.RLeft, epd.SLeft = edgeLeft(fullLeft)
	epd.RHyp, epd.SHyp = edgeHypotenuse(fullHyp)

	epd.CondBottom, epd.CondLeft, epd.CondHyp =
		dfr.checkConditionNumbers(fullBottom, fullLeft, fullHyp)

	return
}

// checkConditionNumbers computes and reports the condition numbers for each edge
func (dfr *DFR2D) checkConditionNumbers(fullBottom, fullLeft,
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

	ABottom := dfr.SolutionBasis.GetInterpMatrix(RmBottom, SmBottom)
	ALeft := dfr.SolutionBasis.GetInterpMatrix(RmLeft, SmLeft)
	AHyp := dfr.SolutionBasis.GetInterpMatrix(RmHyp, SmHyp)

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

// lebesgueConstant computes the Lebesgue constant for a set of interpolation points.
// The Lebesgue constant provides a measure of how good a set of points is for polynomial
// interpolation. It bounds the ratio of the interpolation error to the best approximation error.
func lebesgueConstant(points []float64) float64 {
	n := len(points)
	if n < 2 {
		return 0.0
	}

	// Number of evaluation points
	const numEval = 1000
	maxLebesgue := 0.0

	for i := 0; i < numEval; i++ {
		x := -1.0 + 2.0*float64(i)/float64(numEval-1)
		sum := 0.0

		for j := 0; j < n; j++ {
			// Compute the Lagrange basis polynomial
			basis := 1.0
			for k := 0; k < n; k++ {
				if k != j {
					basis *= (x - points[k]) / (points[j] - points[k])
				}
			}
			sum += math.Abs(basis)
		}

		if sum > maxLebesgue {
			maxLebesgue = sum
		}
	}

	return maxLebesgue
}

func optimizeFunc(y []float64, dfr *DFR2D, edgeFunc func([]float64) (R, S []float64)) float64 {
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
	lebConst := lebesgueConstant(fullX)

	// Return just the Lebesgue constant as our objective
	// We'll separately evaluate and report condition numbers but not use them for optimization
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

// optimizeEdge performs the optimization for a given edge. Since the
// endpoints are fixed, we optimize only the interior points. We transform
// initial interior values from x-space to y-space via atanh.
func optimizeEdge(dfr *DFR2D, initialInterior []float64, edgeFunc func([]float64) (R, S []float64)) []float64 {
	n := len(initialInterior)
	initialY := make([]float64, n)
	for i, xi := range initialInterior {
		// Keep values slightly away from ±1 to avoid numerical issues with atanh
		safeX := math.Max(-0.9999, math.Min(0.9999, xi))
		initialY[i] = math.Atanh(safeX)
	}

	// First, evaluate the condition number of the initial distribution
	fullInitial := make([]float64, n+2)
	fullInitial[0] = -1
	copy(fullInitial[1:n+1], initialInterior)
	fullInitial[n+1] = 1

	// Calculate initial Lebesgue constant
	initLebConst := lebesgueConstant(fullInitial)

	// Calculate initial condition number
	// R, S := edgeFunc(fullInitial)
	// Rm, Sm := utils.NewVector(len(R), R), utils.NewVector(len(S), S)
	// A := dfr.SolutionBasis.GetInterpMatrix(Rm, Sm)
	// initCond := conditionNumber(A.M)

	// fmt.Printf("  Initial distribution - Lebesgue: %.6f, Condition: %.6e\n", initLebConst, initCond)

	problem := optimize.Problem{
		Func: func(y []float64) float64 {
			return optimizeFunc(y, dfr, edgeFunc)
		},
	}

	settings := optimize.Settings{
		MajorIterations:   1000, // Increase max iterations
		GradientThreshold: 1e-6, // Tighter convergence criteria
		Concurrent:        4,    // Use multiple goroutines if possible
	}

	result, err := optimize.Minimize(problem, initialY, &settings, nil)
	if err != nil {
		fmt.Println("  Optimization error:", err)
		// Try a different starting point with LGL points if optimization fails
		if n+2 <= 10 { // Use precalculated LGL points for small n
			lglPoints := legendreGaussLobattoPoints(n + 2)
			// Extract interior points
			lglInterior := make([]float64, n)
			copy(lglInterior, lglPoints[1:n+1])

			// Transform to unconstrained space
			for i, xi := range lglInterior {
				safeX := math.Max(-0.9999, math.Min(0.9999, xi))
				initialY[i] = math.Atanh(safeX)
			}

			// Try optimization again
			result, err = optimize.Minimize(problem, initialY, &settings, nil)

			if err != nil {
				// If still failing, return the LGL points which should be pretty good
				fmt.Println("  Optimization failed. Using LGL points.")
				return lglInterior
			}
		} else {
			// Fall back to initial points
			fmt.Println("  Optimization failed. Using initial points.")
			return initialInterior
		}
	}

	// Transform the optimized y back to x for the interior.
	optXInterior := make([]float64, n)
	for i, yi := range result.X {
		optXInterior[i] = math.Tanh(yi)
	}

	// Check if optimization improved the Lebesgue constant
	fullOptimized := make([]float64, n+2)
	fullOptimized[0] = -1
	copy(fullOptimized[1:n+1], optXInterior)
	fullOptimized[n+1] = 1

	optLebConst := lebesgueConstant(fullOptimized)

	// Calculate optimized condition number
	// R, S = edgeFunc(fullOptimized)
	// Rm, Sm = utils.NewVector(len(R), R), utils.NewVector(len(S), S)
	// A = dfr.SolutionBasis.GetInterpMatrix(Rm, Sm)
	// optCond := conditionNumber(A.M)

	// fmt.Printf("  Optimized distribution - Lebesgue: %.6f, Condition: %.6e\n", optLebConst, optCond)

	// Accept the optimization results even if condition number is high
	// as long as we improved the Lebesgue constant
	if optLebConst > initLebConst*1.01 { // Allow 1% tolerance
		fmt.Printf("  Optimization increased Lebesgue constant (%.4f -> %.4f). Reverting.\n",
			initLebConst, optLebConst)
		return initialInterior
	}

	return optXInterior
}

package DG2D

import (
	"github.com/notargets/gocfd/utils"
	"math"
)

// ModalTransferEdgeInterpolation performs edge interpolation using the modal transfer method
// This is now the recommended approach based on optimization results
type ModalTransferEdgeInterpolation struct {
	le               *LagrangeElement2D
	rt               *RTElement
	transferMatrices [3]utils.Matrix   // Pre-computed transfer matrices for each edge
	edge1DBases      [3]*JacobiBasis1D // 1D bases for each edge
}

// NewModalTransferEdgeInterpolation creates and initializes the modal transfer system
func NewModalTransferEdgeInterpolation(le *LagrangeElement2D, rt *RTElement) *ModalTransferEdgeInterpolation {
	mt := &ModalTransferEdgeInterpolation{
		le: le,
		rt: rt,
	}

	// Pre-compute transfer matrices for each edge
	for edgeNum := 0; edgeNum < 3; edgeNum++ {
		mt.transferMatrices[edgeNum] = mt.buildTransferMatrix(edgeNum)
		mt.edge1DBases[edgeNum] = mt.buildEdge1DBasis(edgeNum)
	}

	return mt
}

// InterpolateToEdge performs the modal transfer interpolation
func (mt *ModalTransferEdgeInterpolation) InterpolateToEdge(edgeNum int,
	interiorValues []float64) ([]float64, error) {

	// Step 1: Convert interior nodal values to 2D modal coefficients
	interiorMatrix := utils.NewMatrix(len(interiorValues), 1, interiorValues)
	modalCoeffs2D := mt.le.JB2D.Vinv.Mul(interiorMatrix)

	// Step 2: Apply modal transfer matrix to get 1D edge modal coefficients
	modalCoeffs1D := mt.transferMatrices[edgeNum].Mul(modalCoeffs2D)

	// Step 3: Convert 1D modal coefficients to nodal values on edge
	edgeNodalMatrix := mt.edge1DBases[edgeNum].V.Mul(modalCoeffs1D)

	// Extract values
	edgeValues := make([]float64, mt.rt.NpEdge)
	for i := 0; i < mt.rt.NpEdge; i++ {
		edgeValues[i] = edgeNodalMatrix.At(i, 0)
	}

	return edgeValues, nil
}

// buildTransferMatrix constructs the modal transfer matrix for a specific edge
func (mt *ModalTransferEdgeInterpolation) buildTransferMatrix(edgeNum int) utils.Matrix {
	num2DModes := mt.le.Np
	num1DModes := mt.rt.NpEdge

	// Create transfer matrix T such that: modal_1D = T * modal_2D
	T := utils.NewMatrix(num1DModes, num2DModes)

	// For each 1D mode
	for m := 0; m < num1DModes; m++ {
		// For each 2D mode
		for j := 0; j < num2DModes; j++ {
			// Get 2D mode orders
			order2D := mt.le.JB2D.Order2DAtJ[j]
			i, k := order2D[0], order2D[1]

			// Compute coupling coefficient
			coeff := mt.computeModalCoupling(i, k, m, edgeNum)
			T.Set(m, j, coeff)
		}
	}

	return T
}

// computeModalCoupling computes the coupling between 2D mode (i,k) and 1D mode m
func (mt *ModalTransferEdgeInterpolation) computeModalCoupling(i, k, m int, edgeNum int) float64 {
	// Use high-order quadrature for accurate integration
	nQuad := 20
	quadPoints, quadWeights := getLegendreQuadrature(nQuad)

	integral := 0.0

	for q := 0; q < nQuad; q++ {
		xi := quadPoints[q]
		weight := quadWeights[q]

		// Map 1D coordinate to 2D edge
		r, s := mapEdgeParamToRS(xi, edgeNum)

		// Evaluate 2D mode at edge point
		val2D := mt.evaluateMode2D(i, k, r, s)

		// Evaluate 1D mode
		val1D := evaluateJacobi1D(xi, m, 0.0, 0.0)

		integral += weight * val2D * val1D
	}

	// Normalize by 1D mode norm
	norm1D := math.Sqrt(2.0 / (2.0*float64(m) + 1.0))

	return integral / norm1D
}

// buildEdge1DBasis creates the 1D Jacobi basis for edge interpolation
func (mt *ModalTransferEdgeInterpolation) buildEdge1DBasis(edgeNum int) *JacobiBasis1D {
	// Get edge points in canonical 1D coordinates
	edgePoints := getEdgePoints(mt.rt, edgeNum)
	edge1DCoords := make([]float64, len(edgePoints))

	for i, pt := range edgePoints {
		edge1DCoords[i] = getCanonicalEdgeParameter(pt, edgeNum)
	}

	// Create 1D Jacobi basis
	edge1DCoordVec := utils.NewVector(len(edge1DCoords), edge1DCoords)
	return NewJacobiBasis1D(mt.rt.NpEdge-1, edge1DCoordVec)
}

// ShockCapturingFilter applies modal filtering for shock capturing
type ShockCapturingFilter struct {
	mt           *ModalTransferEdgeInterpolation
	filterOrders []int   // Modal orders to filter
	filterAlpha  float64 // Filter strength (0 = no filter, 1 = complete removal)
}

// NewShockCapturingFilter creates a filter tuned for shock capturing
func NewShockCapturingFilter(mt *ModalTransferEdgeInterpolation, order int) *ShockCapturingFilter {
	// Filter high-order modes more aggressively
	filterOrders := make([]int, 0)
	cutoff := int(float64(order) * 0.7) // Keep 70% of modes unfiltered

	for i := cutoff; i <= order; i++ {
		filterOrders = append(filterOrders, i)
	}

	return &ShockCapturingFilter{
		mt:           mt,
		filterOrders: filterOrders,
		filterAlpha:  0.0, // Will be set adaptively based on shock indicator
	}
}

// ApplyFilter filters the solution to reduce oscillations near shocks
func (sf *ShockCapturingFilter) ApplyFilter(elementValues []float64,
	shockIndicator float64) []float64 {

	// Adaptive filter strength based on shock indicator
	sf.filterAlpha = math.Min(shockIndicator, 1.0)

	if sf.filterAlpha < 1e-10 {
		return elementValues // No filtering needed
	}

	// Convert to modal space
	nodalMatrix := utils.NewMatrix(len(elementValues), 1, elementValues)
	modalCoeffs := sf.mt.le.JB2D.Vinv.Mul(nodalMatrix)

	// Apply exponential filter to high modes
	for j := 0; j < sf.mt.le.Np; j++ {
		order2D := sf.mt.le.JB2D.Order2DAtJ[j]
		totalOrder := order2D[0] + order2D[1]

		// Check if this mode should be filtered
		for _, filterOrder := range sf.filterOrders {
			if totalOrder >= filterOrder {
				// Exponential filter
				sigma := math.Exp(-sf.filterAlpha * float64(totalOrder-filterOrder+1))
				modalCoeffs.Set(j, 0, modalCoeffs.At(j, 0)*sigma)
				break
			}
		}
	}

	// Convert back to nodal space
	filteredMatrix := sf.mt.le.JB2D.V.Mul(modalCoeffs)

	// Extract filtered values
	filtered := make([]float64, len(elementValues))
	for i := 0; i < len(filtered); i++ {
		filtered[i] = filteredMatrix.At(i, 0)
	}

	return filtered
}

// PerformanceMetrics tracks the performance of modal transfer
type PerformanceMetrics struct {
	ConditionNumber       float64
	MaxInterpolationError float64
	LebesgueConstant      float64
	ComputationTime       float64
}

// EvaluatePerformance computes performance metrics for the modal transfer method
func (mt *ModalTransferEdgeInterpolation) EvaluatePerformance(edgeNum int) PerformanceMetrics {
	metrics := PerformanceMetrics{}

	// Condition number of transfer matrix
	metrics.ConditionNumber = mt.transferMatrices[edgeNum].ConditionNumber()

	// Lebesgue constant for edge points
	edgePoints := getEdgePoints(mt.rt, edgeNum)
	metrics.LebesgueConstant = computeLebesgueConstantStable(edgePoints, edgeNum)

	// Test interpolation accuracy on polynomials
	maxError := 0.0
	testPolys := []func(r, s float64) float64{
		func(r, s float64) float64 { return 1.0 },
		func(r, s float64) float64 { return r },
		func(r, s float64) float64 { return s },
		func(r, s float64) float64 { return r * s },
		func(r, s float64) float64 { return r*r - s*s },
	}

	for _, poly := range testPolys {
		// Evaluate at interior points
		interiorVals := make([]float64, mt.le.Np)
		for i := 0; i < mt.le.Np; i++ {
			r, s := mt.le.R.AtVec(i), mt.le.S.AtVec(i)
			interiorVals[i] = poly(r, s)
		}

		// Interpolate to edge
		edgeVals, _ := mt.InterpolateToEdge(edgeNum, interiorVals)

		// Compute error
		for i, pt := range edgePoints {
			exact := poly(pt.R, pt.S)
			error := math.Abs(edgeVals[i] - exact)
			if error > maxError {
				maxError = error
			}
		}
	}

	metrics.MaxInterpolationError = maxError

	return metrics
}

// Helper functions

func evaluateMode2D(jb2d *JacobiBasis2D, i, k int, r, s float64) float64 {
	RVec := utils.NewVector(1, []float64{r})
	SVec := utils.NewVector(1, []float64{s})
	mode := jb2d.Simplex2DP(RVec, SVec, i, k)
	return mode[0]
}

func (mt *ModalTransferEdgeInterpolation) evaluateMode2D(i, k int, r, s float64) float64 {
	return evaluateMode2D(mt.le.JB2D, i, k, r, s)
}

func evaluateJacobi1D(x float64, n int, alpha, beta float64) float64 {
	// Evaluate Jacobi polynomial P_n^{alpha,beta}(x)
	if n == 0 {
		return 1.0
	}
	if n == 1 {
		return 0.5 * (alpha - beta + (alpha+beta+2.0)*x)
	}

	// Three-term recurrence
	p0 := 1.0
	p1 := 0.5 * (alpha - beta + (alpha+beta+2.0)*x)

	for k := 2; k <= n; k++ {
		a1 := 2.0 * float64(k) * (float64(k) + alpha + beta) * (2.0*float64(k) + alpha + beta - 2.0)
		a2 := (2.0*float64(k) + alpha + beta - 1.0) * (alpha*alpha - beta*beta)
		a3 := (2.0*float64(k) + alpha + beta - 2.0) * (2.0*float64(k) + alpha + beta - 1.0) * (2.0*float64(k) + alpha + beta)
		a4 := 2.0 * (float64(k) + alpha - 1.0) * (float64(k) + beta - 1.0) * (2.0*float64(k) + alpha + beta)

		p2 := ((a2+a3*x)*p1 - a4*p0) / a1
		p0, p1 = p1, p2
	}

	return p1
}

func getLegendreQuadrature(n int) ([]float64, []float64) {
	// Returns Gauss-Legendre quadrature points and weights
	// For production, use a proper quadrature library
	switch n {
	case 20:
		// 20-point Gauss-Legendre quadrature
		points := []float64{
			-0.9931285991850949, -0.9639719272779138, -0.9122344282513259, -0.8391169718222188,
			-0.7463319064601508, -0.6360536807265150, -0.5108670019508271, -0.3737060887154195,
			-0.2277858511416451, -0.0765265211334973, 0.0765265211334973, 0.2277858511416451,
			0.3737060887154195, 0.5108670019508271, 0.6360536807265150, 0.7463319064601508,
			0.8391169718222188, 0.9122344282513259, 0.9639719272779138, 0.9931285991850949,
		}
		weights := []float64{
			0.0176140071391521, 0.0406014298003869, 0.0626720483341091, 0.0832767415767047,
			0.1019301198172404, 0.1181945319615184, 0.1316886384491766, 0.1420961093183820,
			0.1491729864726037, 0.1527533871307258, 0.1527533871307258, 0.1491729864726037,
			0.1420961093183820, 0.1316886384491766, 0.1181945319615184, 0.1019301198172404,
			0.0832767415767047, 0.0626720483341091, 0.0406014298003869, 0.0176140071391521,
		}
		return points, weights
	default:
		// Fallback to lower order
		return getGaussLegendreQuadrature(15)
	}
}

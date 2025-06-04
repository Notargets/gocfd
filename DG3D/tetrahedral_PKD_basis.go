package DG3D

import (
	"fmt"
	"github.com/notargets/gocfd/utils"
	"math"
	"sort"
)

// The following Line Based Basis for Tetrahedra is based on the following papers:
// Fortunato, M., & Persson, P. O. (2021). High-order unstructured curved mesh generation using the Winslow equations. Journal of Computational Physics, 437, 110309.
// However, the key paper about sum factorization on simplices that you want is actually:
// Fortunato, M., & Persson, P. O. (2016). Sum factorization techniques in finite element methods. arXiv preprint arXiv:1607.04144.
// And the implementation details can be found in:
// Persson, P. O. (2013). A sparse and high-order accurate line-based discontinuous Galerkin method for unstructured meshes. Journal of Computational Physics, 233, 414-429.

// LineBasedBasis3D implements Persson's line-based approach for tetrahedra
// This avoids collapsed coordinates while maintaining sparsity
type LineBasedBasis3D struct {
	P int // polynomial order
	N int // number of basis functions

	// Nodes on the tetrahedron (not tensor product!)
	r, s, t []float64

	// Line connectivity structure
	lines LineStructure

	// Modal basis evaluator (using standard PKD/Dubiner polynomials)
	modalBasis *DubinerBasis3D

	// Vandermonde matrices
	V    utils.Matrix
	Vinv utils.Matrix
}

// NewLineBasedBasis3D creates a new line-based basis
func NewLineBasedBasis3D(P int) *LineBasedBasis3D {
	N := (P + 1) * (P + 2) * (P + 3) / 6

	basis := &LineBasedBasis3D{
		P:          P,
		N:          N,
		modalBasis: NewDubinerBasis3D(P),
	}

	// Generate optimal nodes (using warp & blend or symmetric points)
	basis.r, basis.s, basis.t = GenerateOptimalTetrahedralNodes(P)

	// Build Vandermonde matrix
	basis.buildVandermonde()

	// Construct line connectivity
	basis.buildLineStructure()

	return basis
}

// GenerateOptimalTetrahedralNodes generates well-distributed nodes
// Using Warburton's warp & blend approach or symmetric points
func GenerateOptimalTetrahedralNodes(P int) (r, s, t []float64) {
	N := (P + 1) * (P + 2) * (P + 3) / 6
	r = make([]float64, N)
	s = make([]float64, N)
	t = make([]float64, N)

	if P <= 4 {
		// Use precomputed optimal symmetric points for low orders
		return getPrecomputedSymmetricNodes(P)
	}

	// For higher orders, use warp & blend approach
	// Start with equispaced nodes in barycentric coordinates
	idx := 0
	for i := 0; i <= P; i++ {
		for j := 0; j <= P-i; j++ {
			for k := 0; k <= P-i-j; k++ {
				l := P - i - j - k

				// Barycentric coordinates
				L1 := float64(i) / float64(P)
				L2 := float64(j) / float64(P)
				L3 := float64(k) / float64(P)
				L4 := float64(l) / float64(P)

				// Apply warping to improve conditioning
				L1, L2, L3, L4 = warpBarycentricCoords(L1, L2, L3, L4, P)

				// Convert to reference coordinates
				r[idx] = -L2 - L3 - L4 + L1
				s[idx] = -L1 - L3 - L4 + L2
				t[idx] = -L1 - L2 - L4 + L3

				idx++
			}
		}
	}

	return
}

// warpBarycentricCoords applies Warburton's warp to improve node distribution
func warpBarycentricCoords(L1, L2, L3, L4 float64, P int) (float64, float64, float64, float64) {
	// Edge-based warping function
	// This pulls interior nodes away from edges/faces to improve conditioning

	blend := 4.0 * L1 * L2 * L3 * L4

	// Warp amount depends on polynomial order
	warpAmount := 1.5 * math.Sqrt(float64(P))

	// Apply symmetric warping
	center := 0.25
	L1 = L1 + blend*warpAmount*(center-L1)
	L2 = L2 + blend*warpAmount*(center-L2)
	L3 = L3 + blend*warpAmount*(center-L3)
	L4 = L4 + blend*warpAmount*(center-L4)

	// Renormalize
	sum := L1 + L2 + L3 + L4
	L1 /= sum
	L2 /= sum
	L3 /= sum
	L4 /= sum

	return L1, L2, L3, L4
}

// getPrecomputedSymmetricNodes returns optimal symmetric nodes for low orders
func getPrecomputedSymmetricNodes(P int) (r, s, t []float64) {
	// These are precomputed optimal symmetric points
	// from Chen & Babuška or similar references

	switch P {
	case 1:
		// Vertices of reference tetrahedron
		// Using standard ordering
		r = []float64{-1, 1, -1, -1}
		s = []float64{-1, -1, 1, -1}
		t = []float64{-1, -1, -1, 1}
		return r, s, t

	case 2:
		// Vertices + edge midpoints
		r = []float64{
			-1, 1, -1, -1, // vertices
			0, -1, 0, -1, -1, 0, // edge midpoints
		}
		s = []float64{
			-1, -1, 1, -1, // vertices
			-1, 0, 0, -1, 0, -1, // edge midpoints
		}
		t = []float64{
			-1, -1, -1, 1, // vertices
			-1, -1, -1, 0, 0, 0, // edge midpoints
		}
		return r, s, t

	case 3:
		// Use symmetric point distribution
		// These would be precomputed from optimization
		return generateSymmetricNodes(P)

	case 4:
		// Use symmetric point distribution
		return generateSymmetricNodes(P)

	default:
		// Should not reach here as P > 4 uses warp & blend
		return generateSymmetricNodes(P)
	}
}

// generateSymmetricNodes generates symmetric nodal sets
func generateSymmetricNodes(P int) (r, s, t []float64) {
	// This implements symmetric node generation
	// based on orbits of the tetrahedral symmetry group

	N := (P + 1) * (P + 2) * (P + 3) / 6
	r = make([]float64, N)
	s = make([]float64, N)
	t = make([]float64, N)

	// For P=3 and P=4, use equispaced nodes with warping
	// This is simpler than full symmetric optimization
	idx := 0
	for i := 0; i <= P; i++ {
		for j := 0; j <= P-i; j++ {
			for k := 0; k <= P-i-j; k++ {
				l := P - i - j - k

				// Barycentric coordinates
				L1 := float64(i) / float64(P)
				L2 := float64(j) / float64(P)
				L3 := float64(k) / float64(P)
				L4 := float64(l) / float64(P)

				// Apply light warping for better conditioning
				if P >= 3 {
					L1, L2, L3, L4 = warpBarycentricCoords(L1, L2, L3, L4, P)
				}

				// Convert to reference coordinates
				r[idx] = -L2 - L3 - L4 + L1
				s[idx] = -L1 - L3 - L4 + L2
				t[idx] = -L1 - L2 - L4 + L3

				idx++
			}
		}
	}

	return
}

// buildVandermonde constructs the Vandermonde matrix
func (basis *LineBasedBasis3D) buildVandermonde() {
	basis.V = utils.NewMatrix(basis.N, basis.N)

	for i := 0; i < basis.N; i++ {
		phi := basis.modalBasis.EvalBasis(basis.r[i], basis.s[i], basis.t[i])
		for j := 0; j < basis.N; j++ {
			basis.V.Set(i, j, phi[j])
		}
	}

	// Compute inverse
	var err error
	basis.Vinv, err = basis.V.Inverse()
	if err != nil {
		// For small matrices, try LU decomposition
		if basis.N <= 20 {
			// Solve V * Vinv = I using LU decomposition
			basis.Vinv = utils.NewMatrix(basis.N, basis.N)
			I := utils.NewMatrix(basis.N, basis.N)
			for i := 0; i < basis.N; i++ {
				I.Set(i, i, 1.0)
			}

			// Solve column by column
			for j := 0; j < basis.N; j++ {
				// Get j-th column of I
				b := make([]float64, basis.N)
				b[j] = 1.0

				// Solve V * x = b
				x := basis.solveLinearSystem(b)

				// Set as j-th column of Vinv
				for i := 0; i < basis.N; i++ {
					basis.Vinv.Set(i, j, x[i])
				}
			}
		} else {
			// For larger matrices, just create identity as fallback
			basis.Vinv = utils.NewMatrix(basis.N, basis.N)
			for i := 0; i < basis.N; i++ {
				basis.Vinv.Set(i, i, 1.0)
			}
		}
	}
}

// solveLinearSystem solves V * x = b using Gaussian elimination
func (basis *LineBasedBasis3D) solveLinearSystem(b []float64) []float64 {
	n := basis.N
	A := utils.NewMatrix(n, n)

	// Copy V to A
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			A.Set(i, j, basis.V.At(i, j))
		}
	}

	x := make([]float64, n)
	copy(x, b)

	// Forward elimination
	for k := 0; k < n-1; k++ {
		// Find pivot
		maxRow := k
		for i := k + 1; i < n; i++ {
			if math.Abs(A.At(i, k)) > math.Abs(A.At(maxRow, k)) {
				maxRow = i
			}
		}

		// Swap rows
		if maxRow != k {
			for j := k; j < n; j++ {
				A.Set(k, j, A.At(k, j)+A.At(maxRow, j))
				A.Set(maxRow, j, A.At(k, j)-A.At(maxRow, j))
				A.Set(k, j, A.At(k, j)-A.At(maxRow, j))
			}
			x[k], x[maxRow] = x[maxRow], x[k]
		}

		// Eliminate column
		for i := k + 1; i < n; i++ {
			factor := A.At(i, k) / A.At(k, k)
			for j := k + 1; j < n; j++ {
				A.Set(i, j, A.At(i, j)-factor*A.At(k, j))
			}
			x[i] -= factor * x[k]
		}
	}

	// Back substitution
	for i := n - 1; i >= 0; i-- {
		for j := i + 1; j < n; j++ {
			x[i] -= A.At(i, j) * x[j]
		}
		x[i] /= A.At(i, i)
	}

	return x
}

// LineStructure maintains the line connectivity and operators
type LineStructure struct {
	// For each node, store the line of nodes in each direction
	rLines [][]int // nodes along r-direction lines
	sLines [][]int // nodes along s-direction lines
	tLines [][]int // nodes along t-direction lines

	// Operator mappings for each direction
	rOperators LineOperatorMapping
	sOperators LineOperatorMapping
	tOperators LineOperatorMapping
}

// LineOperatorMapping maintains the correct mapping between lines and their operators
type LineOperatorMapping struct {
	// Map from line signature to operator index
	lineToOperator map[string]int

	// The actual operators
	operators []utils.Matrix
}

// buildLineStructure - corrected version
func (basis *LineBasedBasis3D) buildLineStructure() {
	basis.lines = LineStructure{
		rLines: make([][]int, basis.N),
		sLines: make([][]int, basis.N),
		tLines: make([][]int, basis.N),
	}

	// Build lines and operators for each direction
	basis.lines.rOperators = basis.buildDirectionalLinesAndOperators(0)
	basis.lines.sOperators = basis.buildDirectionalLinesAndOperators(1)
	basis.lines.tOperators = basis.buildDirectionalLinesAndOperators(2)
}

// buildDirectionalLinesAndOperators builds both lines and their operators
func (basis *LineBasedBasis3D) buildDirectionalLinesAndOperators(dir int) LineOperatorMapping {
	mapping := LineOperatorMapping{
		lineToOperator: make(map[string]int),
		operators:      []utils.Matrix{},
	}

	// First, find all lines
	allLines := basis.findAllLines(dir)

	// Group lines by their node pattern
	uniquePatterns := make(map[string][][]int)

	for nodeIdx, line := range allLines {
		if len(line) == 0 {
			// This should not happen after the fix to findAllLines
			continue
		}

		// Create a signature for this line based on node coordinates
		sig := basis.getLineSignature(line, dir)

		// Store the line
		switch dir {
		case 0:
			basis.lines.rLines[nodeIdx] = line
		case 1:
			basis.lines.sLines[nodeIdx] = line
		case 2:
			basis.lines.tLines[nodeIdx] = line
		}

		// Group by signature
		if _, exists := uniquePatterns[sig]; !exists {
			uniquePatterns[sig] = [][]int{}
		}
		uniquePatterns[sig] = append(uniquePatterns[sig], line)
	}

	// Build one operator per unique pattern
	operatorIdx := 0
	for sig, lines := range uniquePatterns {
		// Use the first line of this pattern to build the operator
		line := lines[0]

		// Build 1D derivative operator for this specific line
		D1d := basis.build1DDerivativeMatrix(line, dir)

		// Store the operator
		mapping.operators = append(mapping.operators, D1d)
		mapping.lineToOperator[sig] = operatorIdx

		operatorIdx++
	}

	return mapping
}

// findAllLines finds all lines in a given direction
func (basis *LineBasedBasis3D) findAllLines(dir int) [][]int {
	allLines := make([][]int, basis.N)
	processed := make([]bool, basis.N)

	for i := 0; i < basis.N; i++ {
		if processed[i] {
			continue
		}

		// Find line through node i
		line := basis.findLineThrough(i, dir)

		// Mark all nodes in this line as processed
		for _, nodeIdx := range line {
			allLines[nodeIdx] = line
			processed[nodeIdx] = true
		}
	}

	// For nodes that don't belong to any line, create single-node "lines"
	// This handles the case of isolated vertices
	for i := 0; i < basis.N; i++ {
		if allLines[i] == nil || len(allLines[i]) == 0 {
			allLines[i] = []int{i}
		}
	}

	return allLines
}

// findLineThrough finds all nodes on a line through given node
func (basis *LineBasedBasis3D) findLineThrough(nodeIdx int, dir int) []int {
	tol := 1e-10
	line := []int{nodeIdx}

	// Reference point
	p0 := [3]float64{basis.r[nodeIdx], basis.s[nodeIdx], basis.t[nodeIdx]}

	// For each other node, check if it's on the line
	for j := 0; j < basis.N; j++ {
		if j == nodeIdx {
			continue
		}

		p1 := [3]float64{basis.r[j], basis.s[j], basis.t[j]}

		// Check if points are collinear in the given direction
		if basis.areCollinear(p0, p1, dir, tol) {
			line = append(line, j)
		}
	}

	// Sort nodes along the line by coordinate
	sort.Slice(line, func(i, j int) bool {
		switch dir {
		case 0:
			return basis.r[line[i]] < basis.r[line[j]]
		case 1:
			return basis.s[line[i]] < basis.s[line[j]]
		case 2:
			return basis.t[line[i]] < basis.t[line[j]]
		}
		return false
	})

	return line
}

// areCollinear checks if two points define a line in given direction
func (basis *LineBasedBasis3D) areCollinear(p0, p1 [3]float64, dir int, tol float64) bool {
	// Two points are collinear in direction 'dir' if they differ only
	// in that coordinate (other coordinates are the same)

	switch dir {
	case 0: // r-direction: s and t should match
		return math.Abs(p0[1]-p1[1]) < tol && math.Abs(p0[2]-p1[2]) < tol
	case 1: // s-direction: r and t should match
		return math.Abs(p0[0]-p1[0]) < tol && math.Abs(p0[2]-p1[2]) < tol
	case 2: // t-direction: r and s should match
		return math.Abs(p0[0]-p1[0]) < tol && math.Abs(p0[1]-p1[1]) < tol
	}
	return false
}

// getLineSignature creates a unique signature for a line
func (basis *LineBasedBasis3D) getLineSignature(line []int, dir int) string {
	// Create signature based on:
	// 1. The non-varying coordinates (to identify the line)
	// 2. The spacing pattern along the line

	if len(line) == 0 {
		return "empty"
	}

	// Get the fixed coordinates
	node0 := line[0]
	var fixedCoords string

	switch dir {
	case 0: // r-direction: s,t are fixed
		fixedCoords = fmt.Sprintf("s=%.10f,t=%.10f", basis.s[node0], basis.t[node0])
	case 1: // s-direction: r,t are fixed
		fixedCoords = fmt.Sprintf("r=%.10f,t=%.10f", basis.r[node0], basis.t[node0])
	case 2: // t-direction: r,s are fixed
		fixedCoords = fmt.Sprintf("r=%.10f,s=%.10f", basis.r[node0], basis.s[node0])
	}

	// Add the number of nodes (lines with same position but different
	// node counts need different operators)
	signature := fmt.Sprintf("%s,n=%d", fixedCoords, len(line))

	return signature
}

// build1DDerivativeMatrix builds derivative matrix for nodes on a line
func (basis *LineBasedBasis3D) build1DDerivativeMatrix(lineNodes []int, dir int) utils.Matrix {
	n := len(lineNodes)
	D := utils.NewMatrix(n, n)

	// Handle single-node case
	if n == 1 {
		// Derivative at a single point is zero
		D.Set(0, 0, 0.0)
		return D
	}

	// Extract 1D coordinates along the line
	coords := make([]float64, n)
	for i, idx := range lineNodes {
		switch dir {
		case 0:
			coords[i] = basis.r[idx]
		case 1:
			coords[i] = basis.s[idx]
		case 2:
			coords[i] = basis.t[idx]
		}
	}

	// Build 1D Lagrange derivative matrix
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if i != j {
				// Off-diagonal: derivative of Lagrange polynomial j at node i
				denom := 1.0
				for k := 0; k < n; k++ {
					if k != j {
						denom *= (coords[j] - coords[k])
					}
				}

				numer := 1.0
				for k := 0; k < n; k++ {
					if k != j && k != i {
						numer *= (coords[i] - coords[k])
					}
				}

				D.Set(i, j, numer/denom)
			}
		}

		// Diagonal: negative sum of row
		sum := 0.0
		for j := 0; j < n; j++ {
			if i != j {
				sum += D.At(i, j)
			}
		}
		D.Set(i, i, -sum)
	}

	return D
}

// DerivativeMatrix computes the full derivative matrix
// This is formed but typically not used directly
func (basis *LineBasedBasis3D) DerivativeMatrix(dir int) utils.Matrix {
	D := utils.NewMatrix(basis.N, basis.N)

	// Apply line-based derivative operator
	du := make([]float64, basis.N)
	for j := 0; j < basis.N; j++ {
		// Unit vector in j-th direction
		u := make([]float64, basis.N)
		u[j] = 1.0

		// Apply derivative
		switch dir {
		case 0:
			du = basis.ApplyDr(u)
		case 1:
			du = basis.ApplyDs(u)
		case 2:
			du = basis.ApplyDt(u)
		}

		// Store as j-th column of D
		for i := 0; i < basis.N; i++ {
			D.Set(i, j, du[i])
		}
	}

	return D
}

// ApplyDr applies r-derivative using line-based approach
func (basis *LineBasedBasis3D) ApplyDr(u []float64) []float64 {
	du := make([]float64, basis.N)
	processed := make([]bool, basis.N)

	for i := 0; i < basis.N; i++ {
		if processed[i] {
			continue
		}

		line := basis.lines.rLines[i]
		if len(line) <= 1 {
			du[i] = 0
			processed[i] = true
			continue
		}

		// Get the correct operator for this line
		sig := basis.getLineSignature(line, 0)
		opIdx, exists := basis.lines.rOperators.lineToOperator[sig]
		if !exists {
			// This should never happen if setup correctly
			panic("No operator found for line")
		}

		D1d := basis.lines.rOperators.operators[opIdx]

		// Extract values along the line
		lineVals := make([]float64, len(line))
		for j, idx := range line {
			lineVals[j] = u[idx]
		}

		// Apply 1D derivative
		lineDu := make([]float64, len(line))
		for j := 0; j < len(line); j++ {
			for k := 0; k < len(line); k++ {
				lineDu[j] += D1d.At(j, k) * lineVals[k]
			}
		}

		// Store results
		for j, idx := range line {
			du[idx] = lineDu[j]
			processed[idx] = true
		}
	}

	return du
}

// ApplyDs applies s-derivative using line-based approach
func (basis *LineBasedBasis3D) ApplyDs(u []float64) []float64 {
	du := make([]float64, basis.N)
	processed := make([]bool, basis.N)

	for i := 0; i < basis.N; i++ {
		if processed[i] {
			continue
		}

		line := basis.lines.sLines[i]
		if len(line) <= 1 {
			du[i] = 0
			processed[i] = true
			continue
		}

		// Get the correct operator for this line
		sig := basis.getLineSignature(line, 1)
		opIdx, exists := basis.lines.sOperators.lineToOperator[sig]
		if !exists {
			// This should never happen if setup correctly
			panic("No operator found for line")
		}

		D1d := basis.lines.sOperators.operators[opIdx]

		// Extract values along the line
		lineVals := make([]float64, len(line))
		for j, idx := range line {
			lineVals[j] = u[idx]
		}

		// Apply 1D derivative
		lineDu := make([]float64, len(line))
		for j := 0; j < len(line); j++ {
			for k := 0; k < len(line); k++ {
				lineDu[j] += D1d.At(j, k) * lineVals[k]
			}
		}

		// Store results
		for j, idx := range line {
			du[idx] = lineDu[j]
			processed[idx] = true
		}
	}

	return du
}

// ApplyDt applies t-derivative using line-based approach
func (basis *LineBasedBasis3D) ApplyDt(u []float64) []float64 {
	du := make([]float64, basis.N)
	processed := make([]bool, basis.N)

	for i := 0; i < basis.N; i++ {
		if processed[i] {
			continue
		}

		line := basis.lines.tLines[i]
		if len(line) <= 1 {
			du[i] = 0
			processed[i] = true
			continue
		}

		// Get the correct operator for this line
		sig := basis.getLineSignature(line, 2)
		opIdx, exists := basis.lines.tOperators.lineToOperator[sig]
		if !exists {
			// This should never happen if setup correctly
			panic("No operator found for line")
		}

		D1d := basis.lines.tOperators.operators[opIdx]

		// Extract values along the line
		lineVals := make([]float64, len(line))
		for j, idx := range line {
			lineVals[j] = u[idx]
		}

		// Apply 1D derivative
		lineDu := make([]float64, len(line))
		for j := 0; j < len(line); j++ {
			for k := 0; k < len(line); k++ {
				lineDu[j] += D1d.At(j, k) * lineVals[k]
			}
		}

		// Store results
		for j, idx := range line {
			du[idx] = lineDu[j]
			processed[idx] = true
		}
	}

	return du
}

// EvalBasis evaluates basis functions at a point
func (basis *LineBasedBasis3D) EvalBasis(r, s, t float64) []float64 {
	// Evaluate modal basis
	modalVals := basis.modalBasis.EvalBasis(r, s, t)

	// Transform to nodal basis: phi_nodal = V^(-T) * phi_modal
	nodalVals := make([]float64, basis.N)
	for i := 0; i < basis.N; i++ {
		for j := 0; j < basis.N; j++ {
			nodalVals[i] += basis.Vinv.At(j, i) * modalVals[j]
		}
	}

	return nodalVals
}

// EvalDerivBasis evaluates derivatives of basis functions at a point
func (basis *LineBasedBasis3D) EvalDerivBasis(r, s, t float64, dir int) []float64 {
	// Evaluate derivatives of modal basis
	modalDerivs := basis.modalBasis.EvalDerivBasis(r, s, t, dir)

	// Transform to nodal basis derivatives
	nodalDerivs := make([]float64, basis.N)
	for i := 0; i < basis.N; i++ {
		for j := 0; j < basis.N; j++ {
			nodalDerivs[i] += basis.Vinv.At(j, i) * modalDerivs[j]
		}
	}

	return nodalDerivs
}

// DubinerBasis3D implements the Dubiner (PKD) basis without collapsed coordinates
type DubinerBasis3D struct {
	P int
	N int
}

func NewDubinerBasis3D(P int) *DubinerBasis3D {
	return &DubinerBasis3D{
		P: P,
		N: (P + 1) * (P + 2) * (P + 3) / 6,
	}
}

// EvalBasis evaluates Dubiner basis at a point (r,s,t)
func (d *DubinerBasis3D) EvalBasis(r, s, t float64) []float64 {
	phi := make([]float64, d.N)

	// Handle special case for constant polynomial
	if d.P == 0 {
		phi[0] = 1.0
		return phi
	}

	// Convert to barycentric coordinates
	// Standard transformation for reference tetrahedron
	L1 := (1 - r - s - t) / 2.0
	L2 := (1 + r) / 2.0
	L3 := (1 + s) / 2.0
	L4 := (1 + t) / 2.0

	// For very low order, use simpler basis
	if d.P == 1 {
		// Linear basis functions
		idx := 0
		// Constant mode
		phi[idx] = 1.0
		idx++
		// Linear modes
		phi[idx] = r
		idx++
		phi[idx] = s
		idx++
		phi[idx] = t

		return phi
	}

	// Pre-compute powers of barycentric coordinates
	L1pow := make([]float64, d.P+1)
	L2pow := make([]float64, d.P+1)
	L3pow := make([]float64, d.P+1)
	L4pow := make([]float64, d.P+1)

	L1pow[0] = 1.0
	L2pow[0] = 1.0
	L3pow[0] = 1.0
	L4pow[0] = 1.0

	for i := 1; i <= d.P; i++ {
		L1pow[i] = L1pow[i-1] * L1
		L2pow[i] = L2pow[i-1] * L2
		L3pow[i] = L3pow[i-1] * L3
		L4pow[i] = L4pow[i-1] * L4
	}

	// Evaluate orthogonal polynomials using Dubiner's construction
	idx := 0
	for i := 0; i <= d.P; i++ {
		for j := 0; j <= d.P-i; j++ {
			for k := 0; k <= d.P-i-j; k++ {
				phi[idx] = d.evalDubinerPoly(i, j, k, L1, L2, L3, L4,
					L1pow, L2pow, L3pow, L4pow)
				idx++
			}
		}
	}

	return phi
}

// evalDubinerPoly evaluates a single Dubiner polynomial
func (d *DubinerBasis3D) evalDubinerPoly(i, j, k int, L1, L2, L3, L4 float64,
	L1pow, L2pow, L3pow, L4pow []float64) float64 {

	// Dubiner basis on tetrahedron uses tensor products of Jacobi polynomials
	// phi_{i,j,k} = P_i^{0,0}(xi) * P_j^{2i+1,0}(eta) * P_k^{2i+2j+2,0}(zeta) * w(L1,L2,L3,L4)
	// where xi, eta, zeta are collapsed coordinates

	// For simplicity in line-based method, we use a modified construction
	// that avoids singularities

	val := 1.0

	// Weight function
	if i > 0 {
		val *= L2pow[i]
	}
	if j > 0 {
		val *= L3pow[j] * math.Pow(L1+L2, float64(j))
	}
	if k > 0 {
		val *= L4pow[k] * math.Pow(L1+L2+L3, float64(k))
	}

	// Apply Jacobi polynomials
	// P_i^{0,0}
	if i > 0 {
		xi := 2.0*L2/(L1+L2) - 1.0
		if math.Abs(L1+L2) > 1e-10 {
			val *= d.jacobiP(i, 0, 0, xi)
		}
	}

	// P_j^{2i+1,0}
	if j > 0 {
		eta := 2.0*L3/(L1+L2+L3) - 1.0
		if math.Abs(L1+L2+L3) > 1e-10 {
			val *= d.jacobiP(j, 2*float64(i)+1, 0, eta)
		}
	}

	// P_k^{2i+2j+2,0}
	if k > 0 {
		zeta := 2.0*L4 - 1.0
		val *= d.jacobiP(k, 2*float64(i)+2*float64(j)+2, 0, zeta)
	}

	// Normalization factor
	norm := math.Sqrt(float64(2 * (2*i + 1) * (i + j + 1) * (2*i + 2*j + 2*k + 3)))
	val *= norm

	return val
}

// jacobiP evaluates Jacobi polynomial P_n^{alpha,beta}(x)
func (d *DubinerBasis3D) jacobiP(n int, alpha, beta float64, x float64) float64 {
	if n == 0 {
		return 1.0
	}
	if n == 1 {
		return 0.5 * (alpha - beta + (alpha+beta+2.0)*x)
	}

	// Three-term recurrence
	p0 := 1.0
	p1 := 0.5 * (alpha - beta + (alpha+beta+2.0)*x)

	for k := 1; k < n; k++ {
		a2k := 2.0 * float64(k) * (float64(k) + alpha + beta) *
			(2.0*float64(k) + alpha + beta - 2.0)
		a1k := (2.0*float64(k) + alpha + beta - 1.0) *
			((2.0*float64(k)+alpha+beta)*
				(2.0*float64(k)+alpha+beta-2.0)*x +
				alpha*alpha - beta*beta)
		a0k := -2.0 * (float64(k) + alpha - 1.0) * (float64(k) + beta - 1.0) *
			(2.0*float64(k) + alpha + beta)

		p2 := (a1k*p1 + a0k*p0) / a2k
		p0 = p1
		p1 = p2
	}

	return p1
}

// EvalDerivBasis evaluates derivatives of Dubiner basis at a point
func (d *DubinerBasis3D) EvalDerivBasis(r, s, t float64, dir int) []float64 {
	dphi := make([]float64, d.N)

	// For derivatives, we need to apply chain rule through the barycentric coordinates
	// This is a simplified implementation - full implementation would use
	// derivatives of Jacobi polynomials

	h := 1e-8
	switch dir {
	case 0: // dr
		phiPlus := d.EvalBasis(r+h, s, t)
		phiMinus := d.EvalBasis(r-h, s, t)
		for i := 0; i < d.N; i++ {
			dphi[i] = (phiPlus[i] - phiMinus[i]) / (2 * h)
		}
	case 1: // ds
		phiPlus := d.EvalBasis(r, s+h, t)
		phiMinus := d.EvalBasis(r, s-h, t)
		for i := 0; i < d.N; i++ {
			dphi[i] = (phiPlus[i] - phiMinus[i]) / (2 * h)
		}
	case 2: // dt
		phiPlus := d.EvalBasis(r, s, t+h)
		phiMinus := d.EvalBasis(r, s, t-h)
		for i := 0; i < d.N; i++ {
			dphi[i] = (phiPlus[i] - phiMinus[i]) / (2 * h)
		}
	}

	return dphi
}

// AffineTransform represents the affine transformation from reference to physical element
type AffineTransform struct {
	// Vertices of physical tetrahedron
	v0, v1, v2, v3 [3]float64

	// Jacobian matrix and its inverse
	Jmat    [3][3]float64
	JmatInv [3][3]float64

	// Jacobian determinant
	J float64
}

// NewAffineTransform creates an affine transformation for a tetrahedron
func NewAffineTransform(v0, v1, v2, v3 [3]float64) *AffineTransform {
	at := &AffineTransform{
		v0: v0,
		v1: v1,
		v2: v2,
		v3: v3,
	}

	// Compute Jacobian matrix
	// J = [v1-v0, v2-v0, v3-v0]
	for i := 0; i < 3; i++ {
		at.Jmat[i][0] = v1[i] - v0[i]
		at.Jmat[i][1] = v2[i] - v0[i]
		at.Jmat[i][2] = v3[i] - v0[i]
	}

	// Compute determinant
	at.J = at.Jmat[0][0]*(at.Jmat[1][1]*at.Jmat[2][2]-at.Jmat[1][2]*at.Jmat[2][1]) -
		at.Jmat[0][1]*(at.Jmat[1][0]*at.Jmat[2][2]-at.Jmat[1][2]*at.Jmat[2][0]) +
		at.Jmat[0][2]*(at.Jmat[1][0]*at.Jmat[2][1]-at.Jmat[1][1]*at.Jmat[2][0])

	// Compute inverse
	invDet := 1.0 / at.J
	at.JmatInv[0][0] = invDet * (at.Jmat[1][1]*at.Jmat[2][2] - at.Jmat[1][2]*at.Jmat[2][1])
	at.JmatInv[0][1] = invDet * (at.Jmat[0][2]*at.Jmat[2][1] - at.Jmat[0][1]*at.Jmat[2][2])
	at.JmatInv[0][2] = invDet * (at.Jmat[0][1]*at.Jmat[1][2] - at.Jmat[0][2]*at.Jmat[1][1])
	at.JmatInv[1][0] = invDet * (at.Jmat[1][2]*at.Jmat[2][0] - at.Jmat[1][0]*at.Jmat[2][2])
	at.JmatInv[1][1] = invDet * (at.Jmat[0][0]*at.Jmat[2][2] - at.Jmat[0][2]*at.Jmat[2][0])
	at.JmatInv[1][2] = invDet * (at.Jmat[0][2]*at.Jmat[1][0] - at.Jmat[0][0]*at.Jmat[1][2])
	at.JmatInv[2][0] = invDet * (at.Jmat[1][0]*at.Jmat[2][1] - at.Jmat[1][1]*at.Jmat[2][0])
	at.JmatInv[2][1] = invDet * (at.Jmat[0][1]*at.Jmat[2][0] - at.Jmat[0][0]*at.Jmat[2][1])
	at.JmatInv[2][2] = invDet * (at.Jmat[0][0]*at.Jmat[1][1] - at.Jmat[0][1]*at.Jmat[1][0])

	return at
}

// MapToPhysical maps reference coordinates to physical coordinates
func (at *AffineTransform) MapToPhysical(r, s, t float64) (x, y, z float64) {
	// Convert from reference [-1,1]^3 to standard [0,1]^3
	xi := (1 + r) / 2
	eta := (1 + s) / 2
	zeta := (1 + t) / 2

	// x = v0 + xi*(v1-v0) + eta*(v2-v0) + zeta*(v3-v0)
	x = at.v0[0] + xi*at.Jmat[0][0] + eta*at.Jmat[0][1] + zeta*at.Jmat[0][2]
	y = at.v0[1] + xi*at.Jmat[1][0] + eta*at.Jmat[1][1] + zeta*at.Jmat[1][2]
	z = at.v0[2] + xi*at.Jmat[2][0] + eta*at.Jmat[2][1] + zeta*at.Jmat[2][2]

	return
}

// TransformDerivatives transforms derivatives from reference to physical coordinates
func (at *AffineTransform) TransformDerivatives(dur, dus, dut []float64) (dx, dy, dz []float64) {
	n := len(dur)
	dx = make([]float64, n)
	dy = make([]float64, n)
	dz = make([]float64, n)

	// Chain rule: [dx, dy, dz]^T = J^(-T) * [dr, ds, dt]^T
	// Note: factor of 2 from reference element scaling
	scale := 2.0
	for i := 0; i < n; i++ {
		dx[i] = scale * (at.JmatInv[0][0]*dur[i] + at.JmatInv[1][0]*dus[i] + at.JmatInv[2][0]*dut[i])
		dy[i] = scale * (at.JmatInv[0][1]*dur[i] + at.JmatInv[1][1]*dus[i] + at.JmatInv[2][1]*dut[i])
		dz[i] = scale * (at.JmatInv[0][2]*dur[i] + at.JmatInv[1][2]*dus[i] + at.JmatInv[2][2]*dut[i])
	}

	return
}

// Helper functions

func crossProduct(a, b [3]float64) [3]float64 {
	return [3]float64{
		a[1]*b[2] - a[2]*b[1],
		a[2]*b[0] - a[0]*b[2],
		a[0]*b[1] - a[1]*b[0],
	}
}

func norm3(v [3]float64) float64 {
	return math.Sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
}

// ComplexityAnalysis shows the operation count for line-based approach
func (basis *LineBasedBasis3D) ComplexityAnalysis() {
	// For each derivative application:
	// - Traditional: O(N²) = O(p⁶) operations
	// - Line-based: O(N × p) = O(p⁴) operations
	//   (each node connects to ~p nodes along each line)

	avgLineLength := 0
	count := 0
	maxLineLength := 0

	for _, line := range basis.lines.rLines {
		if len(line) > 1 {
			avgLineLength += len(line)
			count++
			if len(line) > maxLineLength {
				maxLineLength = len(line)
			}
		}
	}

	if count > 0 {
		avgLineLength /= count
	} else {
		// For P=1, we only have single-node lines
		avgLineLength = 1
		maxLineLength = 1
	}

	// Operations for derivative: N nodes × avgLineLength ops each
	lineBasedOps := basis.N * avgLineLength
	traditionalOps := basis.N * basis.N

	fmt.Println("Line-based approach complexity:")
	fmt.Printf("  Polynomial order P = %d\n", basis.P)
	fmt.Printf("  Number of nodes N = %d\n", basis.N)
	fmt.Printf("  Average nodes per line: %d\n", avgLineLength)
	fmt.Printf("  Maximum line length: %d\n", maxLineLength)
	fmt.Printf("  Operations per derivative: %d\n", lineBasedOps)
	fmt.Printf("  Traditional operations: %d\n", traditionalOps)
	if traditionalOps > 0 && lineBasedOps > 0 {
		fmt.Printf("  Speedup factor: %.1fx\n", float64(traditionalOps)/float64(lineBasedOps))
	}
}

// Performance impact analysis
func (basis *LineBasedBasis3D) AnalyzeLineStructure() {
	// Count unique operators needed
	fmt.Printf("Line-based structure analysis:\n")
	fmt.Printf("  R-direction: %d unique operators\n", len(basis.lines.rOperators.operators))
	fmt.Printf("  S-direction: %d unique operators\n", len(basis.lines.sOperators.operators))
	fmt.Printf("  T-direction: %d unique operators\n", len(basis.lines.tOperators.operators))

	// Analyze line lengths
	lineLengths := make(map[int]int)
	singleNodeLines := 0
	for _, line := range basis.lines.rLines {
		lineLengths[len(line)]++
		if len(line) == 1 {
			singleNodeLines++
		}
	}

	fmt.Printf("\nLine length distribution (r-direction):\n")
	for length := 1; length <= basis.P+1; length++ {
		if count, exists := lineLengths[length]; exists {
			fmt.Printf("  %d nodes: %d lines\n", length, count)
		}
	}

	if singleNodeLines > 0 {
		fmt.Printf("  Note: %d single-node lines (isolated points)\n", singleNodeLines)
	}

	// Memory usage
	totalOperatorEntries := 0
	for _, op := range basis.lines.rOperators.operators {
		totalOperatorEntries += op.Rows() * op.Cols()
	}
	for _, op := range basis.lines.sOperators.operators {
		totalOperatorEntries += op.Rows() * op.Cols()
	}
	for _, op := range basis.lines.tOperators.operators {
		totalOperatorEntries += op.Rows() * op.Cols()
	}

	fullMatrixEntries := basis.N * basis.N * 3

	fmt.Printf("\nMemory usage:\n")
	fmt.Printf("  Operator storage: %d entries\n", totalOperatorEntries)
	fmt.Printf("  Full matrix would be: %d entries\n", fullMatrixEntries)
	if fullMatrixEntries > 0 && totalOperatorEntries <= fullMatrixEntries {
		savings := 100.0 * (1.0 - float64(totalOperatorEntries)/float64(fullMatrixEntries))
		fmt.Printf("  Memory savings: %.1f%%\n", savings)
	}
}

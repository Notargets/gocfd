package DG3D

import (
	"fmt"
)

// QuadraturePoint represents a quadrature point with coordinates and weight
type QuadraturePoint struct {
	R, S, T float64 // Coordinates in [-1,1]³
	W       float64 // Weight
}

// TetrahedralQuadrature represents a quadrature rule for tetrahedra
type TetrahedralQuadrature struct {
	Order  int
	Points []QuadraturePoint
}

// barycentricToRST converts barycentric coordinates to [-1,1]³ coordinates
// The mapping is based on:
// Unit simplex vertices: (0,0,0), (1,0,0), (0,1,0), (0,0,1)
// Reference tet vertices: (1,-1,-1), (-1,1,-1), (-1,-1,1), (-1,-1,-1)
func barycentricToRST(xi1, xi2, xi3, xi4 float64) (r, s, t float64) {
	r = 2*xi1 - 1
	s = 2*xi2 - 1
	t = 2*xi3 - 1
	return
}

// NewTetrahedralQuadrature creates a Williams-Shunn-Jameson quadrature rule
func NewTetrahedralQuadrature(order int) (*TetrahedralQuadrature, error) {
	switch order {
	case 4:
		return generateOrder4(), nil
	default:
		return nil, fmt.Errorf("order %d not implemented", order)
	}
}

// generateS4 generates a single point at the centroid (all coordinates equal)
func generateS4(value float64) [][4]float64 {
	return [][4]float64{
		{value, value, value, value},
	}
}

// generateS31 generates 4 points with pattern [a, a, a, b]
func generateS31(a, b float64) [][4]float64 {
	return [][4]float64{
		{b, a, a, a},
		{a, b, a, a},
		{a, a, b, a},
		{a, a, a, b},
	}
}

// generateS22 generates 6 points with pattern [a, a, b, b]
func generateS22(a, b float64) [][4]float64 {
	return [][4]float64{
		{a, a, b, b},
		{a, b, a, b},
		{a, b, b, a},
		{b, a, a, b},
		{b, a, b, a},
		{b, b, a, a},
	}
}

// generateS211 generates 12 points with pattern [a, a, b, c]
func generateS211(a, b, c float64) [][4]float64 {
	// For the s211 pattern, we have 12 unique permutations
	// These are all the ways to arrange two a's, one b, and one c
	return [][4]float64{
		{a, a, b, c},
		{a, a, c, b},
		{a, b, a, c},
		{a, b, c, a},
		{a, c, a, b},
		{a, c, b, a},
		{b, a, a, c},
		{b, a, c, a},
		{b, c, a, a},
		{c, a, a, b},
		{c, a, b, a},
		{c, b, a, a},
	}
}

// generateOrder4 creates the 20-point order 4 rule
// Based on Williams-Shunn-Jameson symmetric quadrature rules
func generateOrder4() *TetrahedralQuadrature {
	points := make([]QuadraturePoint, 0, 20)

	// The quadrature rule from Shunn-Ham 2012 for order 4 (degree 4)
	// This integrates polynomials exactly up to degree 4

	// Group 1: 4 points - s31 pattern
	a1 := 0.0323525919372735
	b1 := 0.9029422238811820
	w1 := 0.0070670747944695

	bary1 := generateS31(a1, b1)
	for _, b := range bary1 {
		r, s, t := barycentricToRST(b[0], b[1], b[2], b[3])
		points = append(points, QuadraturePoint{R: r, S: s, T: t, W: w1})
	}

	// Group 2: 12 points - s211 pattern
	a2 := 0.0603604411618854
	b2 := 0.2626825838428973
	c2 := 0.6165965299000219
	w2 := 0.0469986689718877

	bary2 := generateS211(a2, b2, c2)
	for _, b := range bary2 {
		r, s, t := barycentricToRST(b[0], b[1], b[2], b[3])
		points = append(points, QuadraturePoint{R: r, S: s, T: t, W: w2})
	}

	// Group 3: 4 points - s31 pattern
	a3 := 0.3097693042728627
	b3 := 0.0706920867251249
	w3 := 0.0469986689718877

	bary3 := generateS31(a3, b3)
	for _, b := range bary3 {
		r, s, t := barycentricToRST(b[0], b[1], b[2], b[3])
		points = append(points, QuadraturePoint{R: r, S: s, T: t, W: w3})
	}

	// Scale weights to match the [-1,1]³ reference tetrahedron volume
	// The weights from Shunn-Ham are normalized for the unit simplex
	// We need to scale them for our reference element

	// The unit simplex has volume 1/6
	// The [-1,1]³ reference tetrahedron has volume 4/3
	// The transformation Jacobian is 8

	// The weights should sum to the volume of the reference element
	weightSum := 0.0
	for i := range points {
		weightSum += points[i].W
	}

	// Scale factor to make weights sum to 4/3
	scale := ReferenceTetrahedronVolume() / weightSum
	for i := range points {
		points[i].W *= scale
	}

	return &TetrahedralQuadrature{
		Order:  4,
		Points: points,
	}
}

// Integrate performs numerical integration of a function over the reference tetrahedron
func (tq *TetrahedralQuadrature) Integrate(f func(r, s, t float64) float64) float64 {
	sum := 0.0
	for _, pt := range tq.Points {
		sum += pt.W * f(pt.R, pt.S, pt.T)
	}
	return sum
}

// ReferenceTetrahedronVolume returns the volume of the [-1,1]³ reference tetrahedron
func ReferenceTetrahedronVolume() float64 {
	return 4.0 / 3.0
}

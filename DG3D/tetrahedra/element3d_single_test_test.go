package tetrahedra

import (
	"fmt"
	"math"
	"testing"

	"github.com/notargets/gocfd/utils"
)

// TetConfig specifies the configuration for a single tetrahedron
type TetConfig struct {
	Name      string     // Test case name
	Center    [3]float64 // Center point (cx, cy, cz)
	Scale     float64    // Uniform scaling factor
	RotationZ float64    // Rotation angle about Z axis (radians)
	RotationY float64    // Rotation angle about Y axis (radians)
	Shear     float64    // Shear factor (0 for no shear)
}

// buildSingleTetFromConfig creates a single tetrahedron Element3D based on configuration
// This is the primary helper function that all tests should use
func buildSingleTetFromConfig(config TetConfig, N int) *Element3D {
	// Start with reference tetrahedron vertices
	// Reference: (-1,-1,-1), (1,-1,-1), (-1,1,-1), (-1,-1,1)
	refVertices := [][3]float64{
		{-1, -1, -1},
		{1, -1, -1},
		{-1, 1, -1},
		{-1, -1, 1},
	}

	// Apply transformations: scale, shear, rotate, translate
	vertices := make([][3]float64, 4)
	for i := 0; i < 4; i++ {
		// Start with reference coordinates
		x := refVertices[i][0]
		y := refVertices[i][1]
		z := refVertices[i][2]

		// Apply scale
		x *= config.Scale
		y *= config.Scale
		z *= config.Scale

		// Apply shear (in x direction based on y)
		x += config.Shear * y

		// Apply rotation about Z axis
		if config.RotationZ != 0 {
			cosZ := math.Cos(config.RotationZ)
			sinZ := math.Sin(config.RotationZ)
			xNew := cosZ*x - sinZ*y
			yNew := sinZ*x + cosZ*y
			x = xNew
			y = yNew
		}

		// Apply rotation about Y axis
		if config.RotationY != 0 {
			cosY := math.Cos(config.RotationY)
			sinY := math.Sin(config.RotationY)
			xNew := cosY*x + sinY*z
			zNew := -sinY*x + cosY*z
			x = xNew
			z = zNew
		}

		// Apply translation
		vertices[i][0] = x + config.Center[0]
		vertices[i][1] = y + config.Center[1]
		vertices[i][2] = z + config.Center[2]
	}

	// Create single element mesh - VX, VY, VZ are vectors of length 4 (one tet)
	VX := utils.NewVector(4)
	VY := utils.NewVector(4)
	VZ := utils.NewVector(4)
	for i := 0; i < 4; i++ {
		VX.Set(i, vertices[i][0])
		VY.Set(i, vertices[i][1])
		VZ.Set(i, vertices[i][2])
	}

	// Create element-to-vertex connectivity
	EToV := [][]int{{0, 1, 2, 3}}

	// Create TetBasis
	tetBasis := NewTetBasis(N)

	// Initialize Element3D manually since we're creating a test element
	el3d := &Element3D{
		K:        1,
		VX:       VX,
		VY:       VY,
		VZ:       VZ,
		EToV:     EToV,
		TetBasis: tetBasis,
	}

	// Initialize physical coordinates (Np x K matrices)
	el3d.X = utils.NewMatrix(tetBasis.Np, 1)
	el3d.Y = utils.NewMatrix(tetBasis.Np, 1)
	el3d.Z = utils.NewMatrix(tetBasis.Np, 1)

	// Transform reference nodes to physical coordinates
	for n := 0; n < tetBasis.Np; n++ {
		// Get reference coordinates
		r := tetBasis.R.At(n)
		s := tetBasis.S.At(n)
		t := tetBasis.T.At(n)

		// Convert to barycentric coordinates
		L1 := -(1.0 + r + s + t) / 2.0
		L2 := (1.0 + r) / 2.0
		L3 := (1.0 + s) / 2.0
		L4 := (1.0 + t) / 2.0

		// Affine transformation to physical coordinates
		x := L1*vertices[0][0] + L2*vertices[1][0] + L3*vertices[2][0] + L4*vertices[3][0]
		y := L1*vertices[0][1] + L2*vertices[1][1] + L3*vertices[2][1] + L4*vertices[3][1]
		z := L1*vertices[0][2] + L2*vertices[1][2] + L3*vertices[2][2] + L4*vertices[3][2]

		el3d.X.Set(n, 0, x)
		el3d.Y.Set(n, 0, y)
		el3d.Z.Set(n, 0, z)
	}

	// Calculate geometric factors
	el3d.GeometricFactors = el3d.GeometricFactors3D()

	// Calculate face geometry
	el3d.FaceGeometricFactors = el3d.CalcFaceGeometry()

	return el3d
}

// Standard test configurations
func getStandardTestConfigs() []TetConfig {
	return []TetConfig{
		{
			Name:      "Unit at origin",
			Center:    [3]float64{0, 0, 0},
			Scale:     0.5, // Makes unit tet from reference
			RotationZ: 0,
			RotationY: 0,
			Shear:     0,
		},
		{
			Name:      "Translated",
			Center:    [3]float64{2, 3, 4},
			Scale:     0.5,
			RotationZ: 0,
			RotationY: 0,
			Shear:     0,
		},
		{
			Name:      "Scaled",
			Center:    [3]float64{0, 0, 0},
			Scale:     2.0,
			RotationZ: 0,
			RotationY: 0,
			Shear:     0,
		},
		{
			Name:      "Rotated Z=45deg",
			Center:    [3]float64{0, 0, 0},
			Scale:     0.5,
			RotationZ: math.Pi / 4,
			RotationY: 0,
			Shear:     0,
		},
		{
			Name:      "Rotated Y=30deg",
			Center:    [3]float64{0, 0, 0},
			Scale:     0.5,
			RotationZ: 0,
			RotationY: math.Pi / 6,
			Shear:     0,
		},
		{
			Name:      "Double rotation Z=45deg,Y=30deg",
			Center:    [3]float64{0, 0, 0},
			Scale:     0.5,
			RotationZ: math.Pi / 4,
			RotationY: math.Pi / 6,
			Shear:     0,
		},
		{
			Name:      "Sheared",
			Center:    [3]float64{0, 0, 0},
			Scale:     0.5,
			RotationZ: 0,
			RotationY: 0,
			Shear:     0.3,
		},
		{
			Name:      "Complex: translated, scaled, rotated",
			Center:    [3]float64{1, -1, 2},
			Scale:     1.5,
			RotationZ: math.Pi / 3,
			RotationY: math.Pi / 4,
			Shear:     0,
		},
	}
}

// Polynomial orders to test
func getTestPolynomialOrders() []int {
	return []int{1, 2, 3, 4, 5}
}

// Helper functions for geometric computations
func computeTetVolume(v0, v1, v2, v3 [3]float64) float64 {
	// Volume = |det(v1-v0, v2-v0, v3-v0)| / 6
	v10 := [3]float64{v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]}
	v20 := [3]float64{v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]}
	v30 := [3]float64{v3[0] - v0[0], v3[1] - v0[1], v3[2] - v0[2]}

	// Compute determinant
	det := v10[0]*(v20[1]*v30[2]-v20[2]*v30[1]) -
		v10[1]*(v20[0]*v30[2]-v20[2]*v30[0]) +
		v10[2]*(v20[0]*v30[1]-v20[1]*v30[0])

	return math.Abs(det) / 6.0
}

func nearEqual(a, b, tol float64) bool {
	return math.Abs(a-b) < tol
}

// TestSingleTetAffineTransformation validates affine transformations
func TestSingleTetAffineTransformation(t *testing.T) {
	configs := getStandardTestConfigs()
	orders := getTestPolynomialOrders()

	for _, config := range configs {
		for _, N := range orders {
			testName := fmt.Sprintf("%s_N%d", config.Name, N)
			t.Run(testName, func(t *testing.T) {
				el3d := buildSingleTetFromConfig(config, N)

				// Test 1: Check that reference vertices map to expected physical locations
				// Reference vertex coordinates in (r,s,t) space
				refVertexCoords := []struct{ r, s, t float64 }{
					{-1, -1, -1}, // v0
					{1, -1, -1},  // v1
					{-1, 1, -1},  // v2
					{-1, -1, 1},  // v3
				}

				// Find nodes closest to reference vertices
				for i, ref := range refVertexCoords {
					minDist := math.MaxFloat64
					minIdx := -1

					// Find the node closest to this reference vertex
					for j := 0; j < el3d.Np; j++ {
						dr := el3d.R.At(j) - ref.r
						ds := el3d.S.At(j) - ref.s
						dt := el3d.T.At(j) - ref.t
						dist := dr*dr + ds*ds + dt*dt

						if dist < minDist {
							minDist = dist
							minIdx = j
						}
					}

					if minIdx == -1 || minDist > 1e-10 {
						t.Errorf("Vertex %d: no node found at reference position (%.1f,%.1f,%.1f)",
							i, ref.r, ref.s, ref.t)
						continue
					}

					// Check physical coordinates match expected vertex
					expectedX := el3d.VX.At(i)
					expectedY := el3d.VY.At(i)
					expectedZ := el3d.VZ.At(i)

					actualX := el3d.X.At(minIdx, 0)
					actualY := el3d.Y.At(minIdx, 0)
					actualZ := el3d.Z.At(minIdx, 0)

					if !nearEqual(actualX, expectedX, 1e-10) ||
						!nearEqual(actualY, expectedY, 1e-10) ||
						!nearEqual(actualZ, expectedZ, 1e-10) {
						t.Errorf("Vertex %d mapping incorrect: got (%.6f,%.6f,%.6f), want (%.6f,%.6f,%.6f)",
							i, actualX, actualY, actualZ, expectedX, expectedY, expectedZ)
					}
				}

				// Test 2: Validate that the transformation is affine
				// For an affine transformation, the Jacobian should be constant
				J0 := el3d.J.At(0, 0)
				for i := 1; i < el3d.Np; i++ {
					Ji := el3d.J.At(i, 0)
					if !nearEqual(Ji, J0, 1e-10*math.Abs(J0)) {
						t.Errorf("Non-constant Jacobian for affine transformation: J[0]=%f, J[%d]=%f",
							J0, i, Ji)
					}
				}
			})
		}
	}
}

// TestSingleTetJacobian validates Jacobian computations
func TestSingleTetJacobian(t *testing.T) {
	configs := getStandardTestConfigs()
	orders := getTestPolynomialOrders()

	for _, config := range configs {
		for _, N := range orders {
			testName := fmt.Sprintf("%s_N%d", config.Name, N)
			t.Run(testName, func(t *testing.T) {
				el3d := buildSingleTetFromConfig(config, N)

				// Test 1: Jacobian determinant must be positive (valid element)
				for i := 0; i < el3d.Np; i++ {
					J := el3d.J.At(i, 0)
					if J <= 0 {
						t.Errorf("Non-positive Jacobian at node %d: %f", i, J)
					}
				}

				// Test 2: Jacobian relates reference and physical volumes correctly
				// Reference tet with vertices at (-1,-1,-1), (1,-1,-1), (-1,1,-1), (-1,-1,1)
				// has volume = |det(v1-v0, v2-v0, v3-v0)|/6 = |det([2,0,0], [0,2,0], [0,0,2])|/6 = 8/6 = 4/3
				V_ref := 4.0 / 3.0

				// Compute physical volume from vertices
				v0 := [3]float64{el3d.VX.At(0), el3d.VY.At(0), el3d.VZ.At(0)}
				v1 := [3]float64{el3d.VX.At(1), el3d.VY.At(1), el3d.VZ.At(1)}
				v2 := [3]float64{el3d.VX.At(2), el3d.VY.At(2), el3d.VZ.At(2)}
				v3 := [3]float64{el3d.VX.At(3), el3d.VY.At(3), el3d.VZ.At(3)}
				V_phys := computeTetVolume(v0, v1, v2, v3)

				// Expected Jacobian value
				expectedJ := V_phys / V_ref

				// Check Jacobian value (should be constant for affine transformation)
				for i := 0; i < el3d.Np; i++ {
					J := el3d.J.At(i, 0)
					if !nearEqual(J, expectedJ, 1e-10*expectedJ) {
						t.Errorf("Incorrect Jacobian value at node %d: got %f, expected %f",
							i, J, expectedJ)
					}
				}
			})
		}
	}
}

// TestSingleTetMetricTensor validates metric tensor calculations
func TestSingleTetMetricTensor(t *testing.T) {
	configs := getStandardTestConfigs()
	orders := getTestPolynomialOrders()

	for _, config := range configs {
		for _, N := range orders {
			testName := fmt.Sprintf("%s_N%d", config.Name, N)
			t.Run(testName, func(t *testing.T) {
				el3d := buildSingleTetFromConfig(config, N)

				// Test metric tensor identities
				// The metric tensor components (Rx, Ry, Rz, etc.) satisfy:
				// Rx*Xr + Ry*Yr + Rz*Zr = 1
				// Rx*Xs + Ry*Ys + Rz*Zs = 0
				// etc.

				// Compute derivatives of physical coordinates
				Dr, Ds, Dt := el3d.Dr, el3d.Ds, el3d.Dt

				Xr := Dr.Mul(el3d.X)
				Xs := Ds.Mul(el3d.X)
				Xt := Dt.Mul(el3d.X)
				Yr := Dr.Mul(el3d.Y)
				Ys := Ds.Mul(el3d.Y)
				Yt := Dt.Mul(el3d.Y)
				Zr := Dr.Mul(el3d.Z)
				Zs := Ds.Mul(el3d.Z)
				Zt := Dt.Mul(el3d.Z)

				// Check metric tensor identities at each node
				for i := 0; i < el3d.Np; i++ {
					// The metric tensor satisfies: J^(-T) * J^(-1) = I
					// This gives us 9 identities:

					// Row 1: R derivatives
					// Rx*Xr + Ry*Yr + Rz*Zr = 1
					val := el3d.Rx.At(i, 0)*Xr.At(i, 0) +
						el3d.Ry.At(i, 0)*Yr.At(i, 0) +
						el3d.Rz.At(i, 0)*Zr.At(i, 0)
					if !nearEqual(val, 1.0, 1e-10) {
						t.Errorf("Node %d: Rx*Xr + Ry*Yr + Rz*Zr = %f, expected 1.0", i, val)
					}

					// Rx*Xs + Ry*Ys + Rz*Zs = 0
					val = el3d.Rx.At(i, 0)*Xs.At(i, 0) +
						el3d.Ry.At(i, 0)*Ys.At(i, 0) +
						el3d.Rz.At(i, 0)*Zs.At(i, 0)
					if !nearEqual(val, 0.0, 1e-10) {
						t.Errorf("Node %d: Rx*Xs + Ry*Ys + Rz*Zs = %f, expected 0.0", i, val)
					}

					// Rx*Xt + Ry*Yt + Rz*Zt = 0
					val = el3d.Rx.At(i, 0)*Xt.At(i, 0) +
						el3d.Ry.At(i, 0)*Yt.At(i, 0) +
						el3d.Rz.At(i, 0)*Zt.At(i, 0)
					if !nearEqual(val, 0.0, 1e-10) {
						t.Errorf("Node %d: Rx*Xt + Ry*Yt + Rz*Zt = %f, expected 0.0", i, val)
					}

					// Row 2: S derivatives
					// Sx*Xr + Sy*Yr + Sz*Zr = 0
					val = el3d.Sx.At(i, 0)*Xr.At(i, 0) +
						el3d.Sy.At(i, 0)*Yr.At(i, 0) +
						el3d.Sz.At(i, 0)*Zr.At(i, 0)
					if !nearEqual(val, 0.0, 1e-10) {
						t.Errorf("Node %d: Sx*Xr + Sy*Yr + Sz*Zr = %f, expected 0.0", i, val)
					}

					// Sx*Xs + Sy*Ys + Sz*Zs = 1
					val = el3d.Sx.At(i, 0)*Xs.At(i, 0) +
						el3d.Sy.At(i, 0)*Ys.At(i, 0) +
						el3d.Sz.At(i, 0)*Zs.At(i, 0)
					if !nearEqual(val, 1.0, 1e-10) {
						t.Errorf("Node %d: Sx*Xs + Sy*Ys + Sz*Zs = %f, expected 1.0", i, val)
					}

					// Sx*Xt + Sy*Yt + Sz*Zt = 0
					val = el3d.Sx.At(i, 0)*Xt.At(i, 0) +
						el3d.Sy.At(i, 0)*Yt.At(i, 0) +
						el3d.Sz.At(i, 0)*Zt.At(i, 0)
					if !nearEqual(val, 0.0, 1e-10) {
						t.Errorf("Node %d: Sx*Xt + Sy*Yt + Sz*Zt = %f, expected 0.0", i, val)
					}

					// Row 3: T derivatives
					// Tx*Xr + Ty*Yr + Tz*Zr = 0
					val = el3d.Tx.At(i, 0)*Xr.At(i, 0) +
						el3d.Ty.At(i, 0)*Yr.At(i, 0) +
						el3d.Tz.At(i, 0)*Zr.At(i, 0)
					if !nearEqual(val, 0.0, 1e-10) {
						t.Errorf("Node %d: Tx*Xr + Ty*Yr + Tz*Zr = %f, expected 0.0", i, val)
					}

					// Tx*Xs + Ty*Ys + Tz*Zs = 0
					val = el3d.Tx.At(i, 0)*Xs.At(i, 0) +
						el3d.Ty.At(i, 0)*Ys.At(i, 0) +
						el3d.Tz.At(i, 0)*Zs.At(i, 0)
					if !nearEqual(val, 0.0, 1e-10) {
						t.Errorf("Node %d: Tx*Xs + Ty*Ys + Tz*Zs = %f, expected 0.0", i, val)
					}

					// Tx*Xt + Ty*Yt + Tz*Zt = 1
					val = el3d.Tx.At(i, 0)*Xt.At(i, 0) +
						el3d.Ty.At(i, 0)*Yt.At(i, 0) +
						el3d.Tz.At(i, 0)*Zt.At(i, 0)
					if !nearEqual(val, 1.0, 1e-10) {
						t.Errorf("Node %d: Tx*Xt + Ty*Yt + Tz*Zt = %f, expected 1.0", i, val)
					}
				}
			})
		}
	}
}

// TestSingleTetPhysicalDerivatives validates physical derivative transformations
func TestSingleTetPhysicalDerivatives(t *testing.T) {
	configs := getStandardTestConfigs()
	orders := []int{2, 3, 4} // Higher orders for derivative accuracy

	for _, config := range configs {
		for _, N := range orders {
			testName := fmt.Sprintf("%s_N%d", config.Name, N)
			t.Run(testName, func(t *testing.T) {
				el3d := buildSingleTetFromConfig(config, N)

				// Test: Physical derivatives of coordinate functions
				// ∂x/∂x = 1, ∂x/∂y = 0, ∂x/∂z = 0, etc.

				// Compute physical derivatives using the chain rule
				Dr, Ds, Dt := el3d.Dr, el3d.Ds, el3d.Dt

				// ∂u/∂x = Rx*∂u/∂r + Sx*∂u/∂s + Tx*∂u/∂t
				DxX := Dr.Mul(el3d.X).ElMul(el3d.Rx).Add(
					Ds.Mul(el3d.X).ElMul(el3d.Sx)).Add(
					Dt.Mul(el3d.X).ElMul(el3d.Tx))

				DyX := Dr.Mul(el3d.X).ElMul(el3d.Ry).Add(
					Ds.Mul(el3d.X).ElMul(el3d.Sy)).Add(
					Dt.Mul(el3d.X).ElMul(el3d.Ty))

				DzX := Dr.Mul(el3d.X).ElMul(el3d.Rz).Add(
					Ds.Mul(el3d.X).ElMul(el3d.Sz)).Add(
					Dt.Mul(el3d.X).ElMul(el3d.Tz))

				// Check ∂x/∂x ≈ 1, ∂x/∂y ≈ 0, ∂x/∂z ≈ 0
				for i := 0; i < el3d.Np; i++ {
					if !nearEqual(DxX.At(i, 0), 1.0, 1e-8) {
						t.Errorf("∂x/∂x at node %d: got %f, expected 1.0", i, DxX.At(i, 0))
					}
					if !nearEqual(DyX.At(i, 0), 0.0, 1e-8) {
						t.Errorf("∂x/∂y at node %d: got %f, expected 0.0", i, DyX.At(i, 0))
					}
					if !nearEqual(DzX.At(i, 0), 0.0, 1e-8) {
						t.Errorf("∂x/∂z at node %d: got %f, expected 0.0", i, DzX.At(i, 0))
					}
				}

				// Similar tests for Y coordinate
				DxY := Dr.Mul(el3d.Y).ElMul(el3d.Rx).Add(
					Ds.Mul(el3d.Y).ElMul(el3d.Sx)).Add(
					Dt.Mul(el3d.Y).ElMul(el3d.Tx))

				DyY := Dr.Mul(el3d.Y).ElMul(el3d.Ry).Add(
					Ds.Mul(el3d.Y).ElMul(el3d.Sy)).Add(
					Dt.Mul(el3d.Y).ElMul(el3d.Ty))

				DzY := Dr.Mul(el3d.Y).ElMul(el3d.Rz).Add(
					Ds.Mul(el3d.Y).ElMul(el3d.Sz)).Add(
					Dt.Mul(el3d.Y).ElMul(el3d.Tz))

				for i := 0; i < el3d.Np; i++ {
					if !nearEqual(DxY.At(i, 0), 0.0, 1e-8) {
						t.Errorf("∂y/∂x at node %d: got %f, expected 0.0", i, DxY.At(i, 0))
					}
					if !nearEqual(DyY.At(i, 0), 1.0, 1e-8) {
						t.Errorf("∂y/∂y at node %d: got %f, expected 1.0", i, DyY.At(i, 0))
					}
					if !nearEqual(DzY.At(i, 0), 0.0, 1e-8) {
						t.Errorf("∂y/∂z at node %d: got %f, expected 0.0", i, DzY.At(i, 0))
					}
				}
			})
		}
	}
}

// TestSingleTetFaceNormals validates face normal and surface Jacobian calculations
func TestSingleTetFaceNormals(t *testing.T) {
	configs := getStandardTestConfigs()
	orders := getTestPolynomialOrders()

	for _, config := range configs {
		for _, N := range orders {
			testName := fmt.Sprintf("%s_N%d", config.Name, N)
			t.Run(testName, func(t *testing.T) {
				el3d := buildSingleTetFromConfig(config, N)

				// Test face normals
				// For a tetrahedron, we have 4 faces
				// Each face normal should:
				// 1. Point outward
				// 2. Have unit length
				// 3. Be orthogonal to face edges

				Nfaces := 4
				Nfp := el3d.Nfp

				for face := 0; face < Nfaces; face++ {
					// Get nodes on this face
					faceStart := face * Nfp
					faceEnd := (face + 1) * Nfp

					// Check normal consistency across face
					// For affine elements, normals should be constant on each face
					nx0 := el3d.Nx.At(faceStart, 0)
					ny0 := el3d.Ny.At(faceStart, 0)
					nz0 := el3d.Nz.At(faceStart, 0)

					for i := faceStart + 1; i < faceEnd; i++ {
						nx := el3d.Nx.At(i, 0)
						ny := el3d.Ny.At(i, 0)
						nz := el3d.Nz.At(i, 0)

						if !nearEqual(nx, nx0, 1e-10) ||
							!nearEqual(ny, ny0, 1e-10) ||
							!nearEqual(nz, nz0, 1e-10) {
							t.Errorf("Face %d: non-constant normal at node %d", face, i-faceStart)
						}
					}

					// Check that normal is unit length
					// In Hesthaven's implementation, normals are stored as unit vectors
					normalMag := math.Sqrt(nx0*nx0 + ny0*ny0 + nz0*nz0)
					if !nearEqual(normalMag, 1.0, 1e-10) {
						t.Errorf("Face %d: normal not unit length, |n| = %f", face, normalMag)
					}

					// Check surface Jacobian is positive
					for i := faceStart; i < faceEnd; i++ {
						sJ := el3d.SJ.At(i, 0)
						if sJ <= 0 {
							t.Errorf("Face %d: negative surface Jacobian at node %d: %f",
								face, i-faceStart, sJ)
						}
					}

					// Check that surface Jacobian is constant on face (for affine elements)
					sJ0 := el3d.SJ.At(faceStart, 0)
					for i := faceStart + 1; i < faceEnd; i++ {
						sJ := el3d.SJ.At(i, 0)
						if !nearEqual(sJ, sJ0, 1e-10*sJ0) {
							t.Errorf("Face %d: non-constant surface Jacobian", face)
						}
					}
				}
			})
		}
	}
}

// TestSingleTetRobustness tests robustness across different configurations
func TestSingleTetRobustness(t *testing.T) {
	// Test with extreme configurations
	extremeConfigs := []TetConfig{
		{
			Name:      "Very small",
			Center:    [3]float64{0, 0, 0},
			Scale:     0.001,
			RotationZ: 0,
			RotationY: 0,
			Shear:     0,
		},
		{
			Name:      "Very large",
			Center:    [3]float64{0, 0, 0},
			Scale:     1000.0,
			RotationZ: 0,
			RotationY: 0,
			Shear:     0,
		},
		{
			Name:      "High shear",
			Center:    [3]float64{0, 0, 0},
			Scale:     1.0,
			RotationZ: 0,
			RotationY: 0,
			Shear:     0.9,
		},
		{
			Name:      "Near 90 degree rotation",
			Center:    [3]float64{0, 0, 0},
			Scale:     1.0,
			RotationZ: math.Pi/2 - 0.1,
			RotationY: 0,
			Shear:     0,
		},
	}

	for _, config := range extremeConfigs {
		for _, N := range []int{1, 3, 5} {
			testName := fmt.Sprintf("%s_N%d", config.Name, N)
			t.Run(testName, func(t *testing.T) {
				el3d := buildSingleTetFromConfig(config, N)

				// Basic sanity checks
				// 1. Jacobian is positive
				minJ := math.MaxFloat64
				for i := 0; i < el3d.Np; i++ {
					J := el3d.J.At(i, 0)
					if J <= 0 {
						t.Errorf("Non-positive Jacobian at node %d: %f", i, J)
					}
					if J < minJ {
						minJ = J
					}
				}

				// 2. Surface Jacobians are positive
				for i := 0; i < 4*el3d.Nfp; i++ {
					sJ := el3d.SJ.At(i, 0)
					if sJ <= 0 {
						t.Errorf("Non-positive surface Jacobian at face node %d: %f", i, sJ)
					}
				}

				// 3. Check condition number isn't too extreme
				// (This is a rough check - very distorted elements will have large condition numbers)
				if config.Name != "High shear" { // High shear is expected to have poor conditioning
					conditionEstimate := el3d.J.Max() / minJ
					if conditionEstimate > 1e6 {
						t.Logf("Warning: High condition number estimate: %e", conditionEstimate)
					}
				}
			})
		}
	}
}

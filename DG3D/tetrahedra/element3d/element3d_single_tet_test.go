package element3d

import (
	"fmt"
	"math"
	"testing"

	"github.com/notargets/gocfd/DG3D/mesh"
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

// PART 2: Updated buildSingleTetFromConfig that works with the mesh structure
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

	// Create a mesh.Mesh struct with single tetrahedron
	m := mesh.NewMesh()

	// Add vertices with sequential node IDs starting from 1
	for i := 0; i < 4; i++ {
		m.AddNode(i+1, vertices[i][:])
	}

	// Add single tetrahedral element
	// Use node IDs 1-4 (which AddElement will convert to indices 0-3)
	err := m.AddElement(1, mesh.Tet, []int{1}, []int{1, 2, 3, 4})
	if err != nil {
		panic(fmt.Sprintf("Failed to add element: %v", err))
	}

	// Build connectivity
	m.BuildConnectivity()

	// Create Element3D using the new constructor
	el3d, err := NewElement3DFromMesh(N, m)
	if err != nil {
		panic(fmt.Sprintf("Failed to create Element3D from mesh: %v", err))
	}

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

				// Test 1: Check that the transformation maps reference vertices correctly
				// The reference tetrahedron has vertices at:
				refVertices := []struct {
					r, s, t float64
					desc    string
				}{
					{-1, -1, -1, "v0"},
					{1, -1, -1, "v1"},
					{-1, 1, -1, "v2"},
					{-1, -1, 1, "v3"},
				}

				// Compute where these reference vertices should map to in physical space
				// based on the configuration
				expectedPhysicalVertices := make([][3]float64, 4)

				// Start with reference tetrahedron
				refTetVertices := [][3]float64{
					{-1, -1, -1},
					{1, -1, -1},
					{-1, 1, -1},
					{-1, -1, 1},
				}

				// Apply the same transformations we applied in buildSingleTetFromConfig
				for i := 0; i < 4; i++ {
					x := refTetVertices[i][0]
					y := refTetVertices[i][1]
					z := refTetVertices[i][2]

					// Apply scale
					x *= config.Scale
					y *= config.Scale
					z *= config.Scale

					// Apply shear
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
					expectedPhysicalVertices[i][0] = x + config.Center[0]
					expectedPhysicalVertices[i][1] = y + config.Center[1]
					expectedPhysicalVertices[i][2] = z + config.Center[2]
				}

				// Now verify that the Element3D transformation produces the same result
				// For each reference vertex, compute its physical position using barycentric coordinates
				for i, ref := range refVertices {
					// Convert reference coordinates to barycentric
					L1 := -(1.0 + ref.r + ref.s + ref.t) / 2.0
					L2 := (1.0 + ref.r) / 2.0
					L3 := (1.0 + ref.s) / 2.0
					L4 := (1.0 + ref.t) / 2.0

					// The physical position should be the barycentric combination of physical vertices
					// Get the actual physical vertices from VX, VY, VZ
					v0x, v0y, v0z := el3d.VX.At(0), el3d.VY.At(0), el3d.VZ.At(0)
					v1x, v1y, v1z := el3d.VX.At(1), el3d.VY.At(1), el3d.VZ.At(1)
					v2x, v2y, v2z := el3d.VX.At(2), el3d.VY.At(2), el3d.VZ.At(2)
					v3x, v3y, v3z := el3d.VX.At(3), el3d.VY.At(3), el3d.VZ.At(3)

					// Compute the transformed position
					x := L1*v0x + L2*v1x + L3*v2x + L4*v3x
					y := L1*v0y + L2*v1y + L3*v2y + L4*v3y
					z := L1*v0z + L2*v1z + L3*v2z + L4*v3z

					// This should match our expected physical vertex
					if !nearEqual(x, expectedPhysicalVertices[i][0], 1e-10) ||
						!nearEqual(y, expectedPhysicalVertices[i][1], 1e-10) ||
						!nearEqual(z, expectedPhysicalVertices[i][2], 1e-10) {
						t.Errorf("Reference vertex %s maps incorrectly: got (%.6f,%.6f,%.6f), want (%.6f,%.6f,%.6f)",
							ref.desc, x, y, z,
							expectedPhysicalVertices[i][0],
							expectedPhysicalVertices[i][1],
							expectedPhysicalVertices[i][2])
					}

					// Also find the interpolation node closest to this reference vertex
					// and check that it has the correct physical coordinates
					minDist := math.MaxFloat64
					minIdx := -1
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

					// If we found a node at exactly this reference position
					if minDist < 1e-10 {
						actualX := el3d.X.At(minIdx, 0)
						actualY := el3d.Y.At(minIdx, 0)
						actualZ := el3d.Z.At(minIdx, 0)

						// Check it matches the expected physical vertex
						if !nearEqual(actualX, expectedPhysicalVertices[i][0], 1e-10) ||
							!nearEqual(actualY, expectedPhysicalVertices[i][1], 1e-10) ||
							!nearEqual(actualZ, expectedPhysicalVertices[i][2], 1e-10) {
							t.Errorf("Node at reference vertex %s has incorrect physical coordinates: got (%.6f,%.6f,%.6f), want (%.6f,%.6f,%.6f)",
								ref.desc, actualX, actualY, actualZ,
								expectedPhysicalVertices[i][0],
								expectedPhysicalVertices[i][1],
								expectedPhysicalVertices[i][2])
						}
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

				// Test 3: Check that the stored vertices match what we expect
				for i := 0; i < 4; i++ {
					vx := el3d.VX.At(i)
					vy := el3d.VY.At(i)
					vz := el3d.VZ.At(i)

					if !nearEqual(vx, expectedPhysicalVertices[i][0], 1e-10) ||
						!nearEqual(vy, expectedPhysicalVertices[i][1], 1e-10) ||
						!nearEqual(vz, expectedPhysicalVertices[i][2], 1e-10) {
						t.Errorf("Stored vertex %d incorrect: got (%.6f,%.6f,%.6f), want (%.6f,%.6f,%.6f)",
							i, vx, vy, vz,
							expectedPhysicalVertices[i][0],
							expectedPhysicalVertices[i][1],
							expectedPhysicalVertices[i][2])
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

				// Test polynomial exactness
				// For order N element, we should exactly differentiate polynomials up to order N
				t.Run("PolynomialExactness", func(t *testing.T) {
					// Generate test polynomials up to order N
					for p := 0; p <= N; p++ {
						// Generate all monomials of total degree p
						// i.e., all (i,j,k) such that i+j+k = p
						for i := 0; i <= p; i++ {
							for j := 0; j <= p-i; j++ {
								k := p - i - j

								// Skip if we've already tested this monomial
								// (only test a representative subset to avoid too many tests)
								if p > 3 && (i > 0 && j > 0 && k > 0) {
									continue // For high orders, skip some mixed terms
								}

								testName := fmt.Sprintf("x%dy%dz%d", i, j, k)

								// Define the monomial x^i * y^j * z^k
								f := func(x, y, z float64) float64 {
									result := 1.0
									for n := 0; n < i; n++ {
										result *= x
									}
									for n := 0; n < j; n++ {
										result *= y
									}
									for n := 0; n < k; n++ {
										result *= z
									}
									return result
								}

								// Define exact derivatives
								dfdx := func(x, y, z float64) float64 {
									if i == 0 {
										return 0.0
									}
									result := float64(i)
									for n := 0; n < i-1; n++ {
										result *= x
									}
									for n := 0; n < j; n++ {
										result *= y
									}
									for n := 0; n < k; n++ {
										result *= z
									}
									return result
								}

								dfdy := func(x, y, z float64) float64 {
									if j == 0 {
										return 0.0
									}
									result := float64(j)
									for n := 0; n < i; n++ {
										result *= x
									}
									for n := 0; n < j-1; n++ {
										result *= y
									}
									for n := 0; n < k; n++ {
										result *= z
									}
									return result
								}

								dfdz := func(x, y, z float64) float64 {
									if k == 0 {
										return 0.0
									}
									result := float64(k)
									for n := 0; n < i; n++ {
										result *= x
									}
									for n := 0; n < j; n++ {
										result *= y
									}
									for n := 0; n < k-1; n++ {
										result *= z
									}
									return result
								}

								// Evaluate function at nodes
								fVals := utils.NewMatrix(el3d.Np, 1)
								for n := 0; n < el3d.Np; n++ {
									x := el3d.X.At(n, 0)
									y := el3d.Y.At(n, 0)
									z := el3d.Z.At(n, 0)
									fVals.Set(n, 0, f(x, y, z))
								}

								// Compute derivatives
								Dxf := Dr.Mul(fVals).ElMul(el3d.Rx).Add(
									Ds.Mul(fVals).ElMul(el3d.Sx)).Add(
									Dt.Mul(fVals).ElMul(el3d.Tx))

								Dyf := Dr.Mul(fVals).ElMul(el3d.Ry).Add(
									Ds.Mul(fVals).ElMul(el3d.Sy)).Add(
									Dt.Mul(fVals).ElMul(el3d.Ty))

								Dzf := Dr.Mul(fVals).ElMul(el3d.Rz).Add(
									Ds.Mul(fVals).ElMul(el3d.Sz)).Add(
									Dt.Mul(fVals).ElMul(el3d.Tz))

								// Check against exact derivatives
								tol := 1e-10
								if N >= 5 {
									tol = 1e-8 * float64(N) // Scale tolerance with order
								}

								maxErrX, maxErrY, maxErrZ := 0.0, 0.0, 0.0
								for n := 0; n < el3d.Np; n++ {
									x := el3d.X.At(n, 0)
									y := el3d.Y.At(n, 0)
									z := el3d.Z.At(n, 0)

									exactDx := dfdx(x, y, z)
									exactDy := dfdy(x, y, z)
									exactDz := dfdz(x, y, z)

									errX := math.Abs(Dxf.At(n, 0) - exactDx)
									errY := math.Abs(Dyf.At(n, 0) - exactDy)
									errZ := math.Abs(Dzf.At(n, 0) - exactDz)

									maxErrX = math.Max(maxErrX, errX)
									maxErrY = math.Max(maxErrY, errY)
									maxErrZ = math.Max(maxErrZ, errZ)
								}

								// Report errors if they exceed tolerance
								if maxErrX > tol {
									t.Errorf("%s: max ∂f/∂x error = %e (tol = %e)", testName, maxErrX, tol)
								}
								if maxErrY > tol {
									t.Errorf("%s: max ∂f/∂y error = %e (tol = %e)", testName, maxErrY, tol)
								}
								if maxErrZ > tol {
									t.Errorf("%s: max ∂f/∂z error = %e (tol = %e)", testName, maxErrZ, tol)
								}

								if testing.Verbose() && p <= 2 { // Only log low-order results
									t.Logf("%s: errors (∂x=%e, ∂y=%e, ∂z=%e)",
										testName, maxErrX, maxErrY, maxErrZ)
								}
							}
						}
					}
				})
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

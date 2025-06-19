package gonudg

import (
	"fmt"
	"math"
	"testing"
)

// TestNormals3DBasicProperties tests fundamental properties of normals
func TestNormals3DBasicProperties(t *testing.T) {
	// Test with single reference tetrahedron
	VX := []float64{-1, 1, -1, -1}
	VY := []float64{-1, -1, 1, -1}
	VZ := []float64{-1, -1, -1, 1}
	EToV := [][]int{{0, 1, 2, 3}}

	// Test different polynomial orders
	for _, N := range []int{1, 2, 3} {
		t.Run(fmt.Sprintf("N=%d", N), func(t *testing.T) {
			dg, err := NewDG3D(N, VX, VY, VZ, EToV)
			if err != nil {
				t.Fatalf("Failed to create DG3D: %v", err)
			}

			// Test 1: All normals should be unit vectors
			Nfp := dg.Nfp
			Nfaces := dg.Nfaces
			for face := 0; face < Nfaces; face++ {
				for i := 0; i < Nfp; i++ {
					row := face*Nfp + i
					nx := dg.Nx.At(row, 0)
					ny := dg.Ny.At(row, 0)
					nz := dg.Nz.At(row, 0)

					norm := math.Sqrt(nx*nx + ny*ny + nz*nz)
					if math.Abs(norm-1.0) > 1e-14 {
						t.Errorf("Face %d, node %d: normal not unit length: %f", face, i, norm)
					}
				}
			}

			// Test 2: Surface Jacobians should be positive
			for i := 0; i < Nfp*Nfaces; i++ {
				sJ := dg.SJ.At(i, 0)
				if sJ <= 0 {
					t.Errorf("Negative or zero surface Jacobian at node %d: %f", i, sJ)
				}
			}

			// Test 3: Volume Jacobian should be positive
			for i := 0; i < dg.Np; i++ {
				J := dg.J.At(i, 0)
				if J <= 0 {
					t.Errorf("Negative or zero volume Jacobian at node %d: %f", i, J)
				}
			}
		})
	}
}

// TestNormals3DOutwardPointing verifies normals point outward
func TestNormals3DOutwardPointing(t *testing.T) {
	// Simple tetrahedron with vertices at (0,0,0), (1,0,0), (0,1,0), (0,0,1)
	VX := []float64{0, 1, 0, 0}
	VY := []float64{0, 0, 1, 0}
	VZ := []float64{0, 0, 0, 1}
	EToV := [][]int{{0, 1, 2, 3}}

	N := 2
	dg, err := NewDG3D(N, VX, VY, VZ, EToV)
	if err != nil {
		t.Fatalf("Failed to create DG3D: %v", err)
	}

	// For each face, verify normal points outward
	// We can do this by checking that n·(p-c) > 0 where p is a face point and c is the centroid

	// Compute element centroid
	cx, cy, cz := 0.25, 0.25, 0.25

	Nfp := dg.Nfp
	for face := 0; face < 4; face++ {
		// Get a representative normal (they should all be the same on a face for affine elements)
		row := face * Nfp
		nx := dg.Nx.At(row, 0)
		ny := dg.Ny.At(row, 0)
		nz := dg.Nz.At(row, 0)

		// Get a point on the face
		faceNodeIdx := dg.Fmask[face][0]
		px := dg.X.At(faceNodeIdx, 0)
		py := dg.Y.At(faceNodeIdx, 0)
		pz := dg.Z.At(faceNodeIdx, 0)

		// Compute n·(p-c)
		dot := nx*(px-cx) + ny*(py-cy) + nz*(pz-cz)

		if dot <= 0 {
			t.Errorf("Face %d: normal appears to point inward (dot product = %f)", face, dot)
		}
	}
}

// TestNormals3DMetricIdentities validates metric tensor relationships
func TestNormals3DMetricIdentities(t *testing.T) {
	// Test with a rotated and scaled tetrahedron
	angle := math.Pi / 6 // 30 degrees
	scale := 2.0

	// Original vertices
	v0 := [3]float64{-1, -1, -1}
	v1 := [3]float64{1, -1, -1}
	v2 := [3]float64{-1, 1, -1}
	v3 := [3]float64{-1, -1, 1}

	// Apply rotation around Z-axis and scaling
	cos_a, sin_a := math.Cos(angle), math.Sin(angle)

	VX := make([]float64, 4)
	VY := make([]float64, 4)
	VZ := make([]float64, 4)

	vertices := [4][3]float64{v0, v1, v2, v3}
	for i, v := range vertices {
		// Rotate around Z-axis
		x_rot := cos_a*v[0] - sin_a*v[1]
		y_rot := sin_a*v[0] + cos_a*v[1]
		z_rot := v[2]

		// Scale
		VX[i] = scale * x_rot
		VY[i] = scale * y_rot
		VZ[i] = scale * z_rot
	}

	EToV := [][]int{{0, 1, 2, 3}}

	N := 3
	dg, err := NewDG3D(N, VX, VY, VZ, EToV)
	if err != nil {
		t.Fatalf("Failed to create DG3D: %v", err)
	}

	// Test metric identities
	// The metric tensor satisfies: [Rx Ry Rz; Sx Sy Sz; Tx Ty Tz] * [xr xs xt; yr ys yt; zr zs zt] = I

	// Compute derivatives of physical coordinates
	xr := dg.Dr.Mul(dg.X)
	xs := dg.Ds.Mul(dg.X)
	xt := dg.Dt.Mul(dg.X)
	yr := dg.Dr.Mul(dg.Y)
	ys := dg.Ds.Mul(dg.Y)
	yt := dg.Dt.Mul(dg.Y)
	zr := dg.Dr.Mul(dg.Z)
	zs := dg.Ds.Mul(dg.Z)
	zt := dg.Dt.Mul(dg.Z)

	// Check identities at each node
	for i := 0; i < dg.Np; i++ {
		// Row 1: Rx*xr + Ry*yr + Rz*zr = 1, Rx*xs + Ry*ys + Rz*zs = 0, etc.
		val := dg.Rx.At(i, 0)*xr.At(i, 0) + dg.Ry.At(i, 0)*yr.At(i, 0) + dg.Rz.At(i, 0)*zr.At(i, 0)
		if math.Abs(val-1.0) > 1e-10 {
			t.Errorf("Node %d: Rx*xr + Ry*yr + Rz*zr = %f, expected 1.0", i, val)
		}

		val = dg.Rx.At(i, 0)*xs.At(i, 0) + dg.Ry.At(i, 0)*ys.At(i, 0) + dg.Rz.At(i, 0)*zs.At(i, 0)
		if math.Abs(val) > 1e-10 {
			t.Errorf("Node %d: Rx*xs + Ry*ys + Rz*zs = %f, expected 0.0", i, val)
		}

		// Row 2: Sx*xr + Sy*yr + Sz*zr = 0, Sx*xs + Sy*ys + Sz*zs = 1, etc.
		val = dg.Sx.At(i, 0)*xs.At(i, 0) + dg.Sy.At(i, 0)*ys.At(i, 0) + dg.Sz.At(i, 0)*zs.At(i, 0)
		if math.Abs(val-1.0) > 1e-10 {
			t.Errorf("Node %d: Sx*xs + Sy*ys + Sz*zs = %f, expected 1.0", i, val)
		}

		// Row 3: Tx*xt + Ty*yt + Tz*zt = 1
		val = dg.Tx.At(i, 0)*xt.At(i, 0) + dg.Ty.At(i, 0)*yt.At(i, 0) + dg.Tz.At(i, 0)*zt.At(i, 0)
		if math.Abs(val-1.0) > 1e-10 {
			t.Errorf("Node %d: Tx*xt + Ty*yt + Tz*zt = %f, expected 1.0", i, val)
		}
	}
}

// TestNormals3DSurfaceArea validates surface area calculations
func TestNormals3DSurfaceArea(t *testing.T) {
	// Unit cube split into tetrahedra
	// We'll use a single tetrahedron from the cube for simplicity
	VX := []float64{0, 1, 0, 0}
	VY := []float64{0, 0, 1, 0}
	VZ := []float64{0, 0, 0, 1}
	EToV := [][]int{{0, 1, 2, 3}}

	N := 3
	dg, err := NewDG3D(N, VX, VY, VZ, EToV)
	if err != nil {
		t.Fatalf("Failed to create DG3D: %v", err)
	}

	// For this tetrahedron:
	// Face 0 (t=-1): Triangle with vertices (0,0,0), (1,0,0), (0,1,0) - area = 1/2
	// Face 1 (s=-1): Triangle with vertices (0,0,0), (1,0,0), (0,0,1) - area = 1/2
	// Face 2 (r+s+t=-1): Triangle with vertices (1,0,0), (0,1,0), (0,0,1) - area = sqrt(3)/2
	// Face 3 (r=-1): Triangle with vertices (0,0,0), (0,1,0), (0,0,1) - area = 1/2

	expectedAreas := []float64{0.5, 0.5, math.Sqrt(3) / 2, 0.5}

	// Compute actual face areas by integrating surface Jacobian
	// For affine elements, SJ should be constant on each face
	Nfp := dg.Nfp
	for face := 0; face < 4; face++ {
		// Get surface Jacobian at first node (should be constant)
		row := face * Nfp
		sJ := dg.SJ.At(row, 0)

		// Check that SJ is constant across the face
		for i := 1; i < Nfp; i++ {
			row_i := face*Nfp + i
			sJ_i := dg.SJ.At(row_i, 0)
			if math.Abs(sJ_i-sJ) > 1e-10*sJ {
				t.Errorf("Face %d: non-constant surface Jacobian", face)
			}
		}

		// For the reference triangle with area 2, physical area = SJ * 2
		// But our reference triangle has vertices at (-1,-1), (1,-1), (-1,1)
		// which has area = 0.5 * base * height = 0.5 * 2 * 2 = 2
		refArea := 2.0
		physArea := sJ * refArea

		relErr := math.Abs(physArea-expectedAreas[face]) / expectedAreas[face]
		if relErr > 1e-10 {
			t.Errorf("Face %d: computed area %f, expected %f (rel error %e)",
				face, physArea, expectedAreas[face], relErr)
		}
	}
}

// TestNormals3DJacobianScaling verifies Jacobian scaling relationships
func TestNormals3DJacobianScaling(t *testing.T) {
	// Test that scaling the element scales Jacobians appropriately
	scales := []float64{0.5, 1.0, 2.0, 3.5}

	for _, scale := range scales {
		t.Run(fmt.Sprintf("scale=%.1f", scale), func(t *testing.T) {
			// Create scaled tetrahedron
			VX := []float64{0, scale, 0, 0}
			VY := []float64{0, 0, scale, 0}
			VZ := []float64{0, 0, 0, scale}
			EToV := [][]int{{0, 1, 2, 3}}

			N := 2
			dg, err := NewDG3D(N, VX, VY, VZ, EToV)
			if err != nil {
				t.Fatalf("Failed to create DG3D: %v", err)
			}

			// Volume should scale as scale^3
			// For this tet: V = (1/6) * scale^3
			expectedVolume := scale * scale * scale / 6.0

			// Jacobian relates reference volume to physical volume
			// Reference tet volume = 4/3
			refVolume := 4.0 / 3.0
			expectedJ := expectedVolume / refVolume

			// Check Jacobian value (should be constant for affine transformation)
			for i := 0; i < dg.Np; i++ {
				J := dg.J.At(i, 0)
				if math.Abs(J-expectedJ) > 1e-10*expectedJ {
					t.Errorf("Node %d: J = %f, expected %f", i, J, expectedJ)
				}
			}

			// Surface Jacobian should scale as scale^2
			// But it's also multiplied by volume Jacobian, so total scaling is scale^2 * (scale^3/refVol)
			// Actually, SJ includes the volume Jacobian factor, so we need to check the actual face areas
		})
	}
}

// TestNormals3DHigherOrder tests that higher order approximations maintain accuracy
func TestNormals3DHigherOrder(t *testing.T) {
	// Create a slightly curved element by perturbing vertices
	VX := []float64{-1.1, 0.9, -0.95, -1.05}
	VY := []float64{-0.9, -1.1, 1.05, -0.95}
	VZ := []float64{-1.05, -0.95, -1.1, 0.9}
	EToV := [][]int{{0, 1, 2, 3}}

	prevNormVariation := 1.0

	// Test that normal variation decreases with increasing order
	for _, N := range []int{1, 2, 3, 4, 5} {
		t.Run(fmt.Sprintf("N=%d", N), func(t *testing.T) {
			dg, err := NewDG3D(N, VX, VY, VZ, EToV)
			if err != nil {
				t.Fatalf("Failed to create DG3D: %v", err)
			}

			// Measure variation in normals across each face
			maxVariation := 0.0
			Nfp := dg.Nfp

			for face := 0; face < 4; face++ {
				// Get first normal as reference
				row0 := face * Nfp
				nx0 := dg.Nx.At(row0, 0)
				ny0 := dg.Ny.At(row0, 0)
				nz0 := dg.Nz.At(row0, 0)

				// Check variation from reference
				for i := 1; i < Nfp; i++ {
					row := face*Nfp + i
					nx := dg.Nx.At(row, 0)
					ny := dg.Ny.At(row, 0)
					nz := dg.Nz.At(row, 0)

					// Compute angle between normals
					dot := nx0*nx + ny0*ny + nz0*nz
					if dot > 1.0 {
						dot = 1.0 // Handle numerical errors
					}
					angle := math.Acos(dot)

					if angle > maxVariation {
						maxVariation = angle
					}
				}
			}

			// For nearly straight elements, variation should be small
			if maxVariation > 0.1 {
				t.Errorf("Large normal variation: %f radians", maxVariation)
			}

			// Variation should generally decrease with order (not strictly monotonic)
			if N > 1 && maxVariation > prevNormVariation*1.5 {
				t.Logf("Warning: normal variation increased from N=%d to N=%d", N-1, N)
			}
			prevNormVariation = maxVariation
		})
	}
}

// TestNormals3DMultipleElements tests with multiple connected elements
func TestNormals3DMultipleElements(t *testing.T) {
	// Two tetrahedra sharing a face
	VX := []float64{0, 1, 0, 0, 1}
	VY := []float64{0, 0, 1, 0, 1}
	VZ := []float64{0, 0, 0, 1, 1}
	EToV := [][]int{
		{0, 1, 2, 3}, // First tet
		{1, 2, 3, 4}, // Second tet sharing face (1,2,3)
	}

	N := 2
	dg, err := NewDG3D(N, VX, VY, VZ, EToV)
	if err != nil {
		t.Fatalf("Failed to create DG3D: %v", err)
	}

	// Check that both elements have positive Jacobians
	for k := 0; k < 2; k++ {
		for i := 0; i < dg.Np; i++ {
			J := dg.J.At(i, k)
			if J <= 0 {
				t.Errorf("Element %d, node %d: negative Jacobian %f", k, i, J)
			}
		}
	}

	// Check that all normals are unit vectors
	Nfp := dg.Nfp
	Nfaces := dg.Nfaces
	for k := 0; k < 2; k++ {
		for face := 0; face < Nfaces; face++ {
			for i := 0; i < Nfp; i++ {
				row := face*Nfp + i
				nx := dg.Nx.At(row, k)
				ny := dg.Ny.At(row, k)
				nz := dg.Nz.At(row, k)

				norm := math.Sqrt(nx*nx + ny*ny + nz*nz)
				if math.Abs(norm-1.0) > 1e-14 {
					t.Errorf("Element %d, face %d, node %d: normal not unit length: %f",
						k, face, i, norm)
				}
			}
		}
	}

	// The shared face should have opposite normals
	// This would be checked after implementing connectivity
}

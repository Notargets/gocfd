package gonudg

import (
	"fmt"
	"github.com/notargets/gocfd/utils"
	"math"
	"testing"
)

func TestStartUp3D(t *testing.T) {
	// Simple test mesh - single tetrahedron
	VX := []float64{0, 1, 0, 0}
	VY := []float64{0, 0, 1, 0}
	VZ := []float64{0, 0, 0, 1}
	EToV := [][]int{{0, 1, 2, 3}} // 0-based indices

	// Test different polynomial orders
	orders := []int{1, 2, 3, 4}

	for _, N := range orders {
		t.Run(fmt.Sprintf("N=%d", N), func(t *testing.T) {
			dg, err := NewDG3D(N, VX, VY, VZ, EToV)
			if err != nil {
				t.Fatalf("Failed to create DG3D: %v", err)
			}

			// Check basic properties
			expectedNp := (N + 1) * (N + 2) * (N + 3) / 6
			if dg.Np != expectedNp {
				t.Errorf("Wrong Np: got %d, want %d", dg.Np, expectedNp)
			}

			expectedNfp := (N + 1) * (N + 2) / 2
			if dg.Nfp != expectedNfp {
				t.Errorf("Wrong Nfp: got %d, want %d", dg.Nfp, expectedNfp)
			}

			// Check matrix dimensions
			checkMatrixDims := func(name string, mat utils.Matrix, expectedRows, expectedCols int) {
				nr, nc := mat.Dims()
				if nr != expectedRows || nc != expectedCols {
					t.Errorf("%s dimensions wrong: got (%d,%d), want (%d,%d)",
						name, nr, nc, expectedRows, expectedCols)
				}
			}

			checkMatrixDims("V", dg.V, dg.Np, dg.Np)
			checkMatrixDims("Vinv", dg.Vinv, dg.Np, dg.Np)
			checkMatrixDims("MassMatrix", dg.MassMatrix, dg.Np, dg.Np)
			checkMatrixDims("Dr", dg.Dr, dg.Np, dg.Np)
			checkMatrixDims("Ds", dg.Ds, dg.Np, dg.Np)
			checkMatrixDims("Dt", dg.Dt, dg.Np, dg.Np)
			checkMatrixDims("X", dg.X, dg.Np, dg.K)
			checkMatrixDims("Y", dg.Y, dg.Np, dg.K)
			checkMatrixDims("Z", dg.Z, dg.Np, dg.K)

			// Check V*Vinv = I
			I := dg.V.Mul(dg.Vinv)
			nr, nc := I.Dims()
			for i := 0; i < nr; i++ {
				for j := 0; j < nc; j++ {
					expected := 0.0
					if i == j {
						expected = 1.0
					}
					if math.Abs(I.At(i, j)-expected) > 1e-10 {
						t.Errorf("V*Vinv not identity at (%d,%d): got %v",
							i, j, I.At(i, j))
					}
				}
			}

			// Check that physical coordinates are reasonable
			// All nodes should be within the tetrahedron bounds
			for i := 0; i < dg.Np; i++ {
				x := dg.X.At(i, 0)
				y := dg.Y.At(i, 0)
				z := dg.Z.At(i, 0)

				// Check bounds (all coordinates should be >= 0 and X+Y+Z <= 1)
				if x < -1e-10 || y < -1e-10 || z < -1e-10 {
					t.Errorf("Node %d has negative coordinate: (%v,%v,%v)",
						i, x, y, z)
				}
				if x+y+z > 1.0+1e-10 {
					t.Errorf("Node %d outside tetrahedron: X+Y+Z = %v",
						i, x+y+z)
				}
			}

			// Check face masks
			if len(dg.Fmask) != 4 {
				t.Errorf("Wrong number of faces: got %d, want 4", len(dg.Fmask))
			}

			// Each face should have Nfp nodes
			for face := 0; face < 4; face++ {
				if len(dg.Fmask[face]) != dg.Nfp {
					t.Errorf("Face %d has wrong number of nodes: got %d, want %d",
						face, len(dg.Fmask[face]), dg.Nfp)
				}
			}
		})
	}
}

func TestDG3DSimpleMesh(t *testing.T) {
	// Test with a mesh of 2 tetrahedra
	VX := []float64{0, 1, 0, 0, 1}
	VY := []float64{0, 0, 1, 0, 1}
	VZ := []float64{0, 0, 0, 1, 1}
	EToV := [][]int{
		{0, 1, 2, 3}, // First tet
		{1, 2, 3, 4}, // Second tet sharing a face
	}

	N := 2
	dg, err := NewDG3D(N, VX, VY, VZ, EToV)
	if err != nil {
		t.Fatalf("Failed to create DG3D: %v", err)
	}

	// Check that we have 2 elements
	if dg.K != 2 {
		t.Errorf("Wrong number of elements: got %d, want 2", dg.K)
	}

	// Check coordinate matrix dimensions
	nr, nc := dg.X.Dims()
	if nc != 2 {
		t.Errorf("Coordinate matrices should have 2 columns (elements): got %d", nc)
	}
	if nr != dg.Np {
		t.Errorf("Coordinate matrices should have %d rows: got %d", dg.Np, nr)
	}
}

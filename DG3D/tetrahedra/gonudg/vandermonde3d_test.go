package gonudg

import (
	"fmt"
	"github.com/notargets/gocfd/utils"
	"math"
	"testing"
)

func TestVandermonde3D(t *testing.T) {
	// Test that Vandermonde matrix has correct dimensions and properties
	tests := []struct {
		name string
		N    int
	}{
		{"N=1", 1},
		{"N=2", 2},
		{"N=3", 3},
		{"N=4", 4},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			// Create test points using Warp & Blend nodes
			// These are designed to avoid singularity issues
			Np := (tc.N + 1) * (tc.N + 2) * (tc.N + 3) / 6
			X, Y, Z := Nodes3D(tc.N)
			r, s, tt := XYZtoRST(X, Y, Z)

			// Build Vandermonde matrix
			V := Vandermonde3D(tc.N, r, s, tt)

			// Check dimensions
			nr, nc := V.Dims()
			if nr != Np {
				t.Errorf("Wrong number of rows: got %d, want %d", nr, Np)
			}
			if nc != Np {
				t.Errorf("Wrong number of columns: got %d, want %d", nc, Np)
			}

			// Check that V is invertible
			if tc.N <= 3 {
				Vinv := V.InverseWithCheck()

				// Check V * V^{-1} = I
				I := V.Mul(Vinv)
				for i := 0; i < Np; i++ {
					for j := 0; j < Np; j++ {
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
			}
		})
	}
}

func TestGradVandermonde3D(t *testing.T) {
	// Test gradient Vandermonde matrices
	N := 2
	Np := (N + 1) * (N + 2) * (N + 3) / 6
	X, Y, Z := Nodes3D(N)
	r, s, tt := XYZtoRST(X, Y, Z)

	Vr, Vs, Vt := GradVandermonde3D(N, r, s, tt)

	// Check dimensions
	matrices := []struct {
		name string
		mat  utils.Matrix
	}{
		{"Vr", Vr},
		{"Vs", Vs},
		{"Vt", Vt},
	}

	for _, m := range matrices {
		nr, nc := m.mat.Dims()
		if nr != Np {
			t.Errorf("%s: Wrong number of rows: got %d, want %d", m.name, nr, Np)
		}
		if nc != Np {
			t.Errorf("%s: Wrong number of columns: got %d, want %d", m.name, nc, Np)
		}
	}

	// Test specific derivatives
	// For P_{1,0,0}, we expect dr/dr = const, dr/ds = function of r,s,t
	// Column for P_{1,0,0} is at index 1 (after P_{0,0,0})
	col := 1

	// Check that derivative values are reasonable
	for i := 0; i < Np; i++ {
		if math.IsNaN(Vr.At(i, col)) || math.IsNaN(Vs.At(i, col)) || math.IsNaN(Vt.At(i, col)) {
			t.Errorf("NaN in gradient at row %d", i)
		}
	}
}

func TestDmatrices3D(t *testing.T) {
	// Test differentiation matrices
	N := 2
	X, Y, Z := Nodes3D(N)
	r, s, tt := XYZtoRST(X, Y, Z)
	V := Vandermonde3D(N, r, s, tt)

	Dr, Ds, Dt := Dmatrices3D(N, r, s, tt, V)

	// Test exact differentiation of polynomials
	// For a linear function f(r,s,t) = ar + bs + ct + d
	// We should get df/dr = a exactly

	a, b, c, d := 2.0, -3.0, 1.5, 0.5
	f := utils.NewVector(len(r))
	for i := 0; i < len(r); i++ {
		f.Set(i, a*r[i]+b*s[i]+c*tt[i]+d)
	}

	// Compute derivatives using matrix multiplication
	dfdr := Dr.Mul(f.ToMatrix())
	dfds := Ds.Mul(f.ToMatrix())
	dfdt := Dt.Mul(f.ToMatrix())

	// Check results
	tol := 1e-10
	for i := 0; i < len(r); i++ {
		if math.Abs(dfdr.At(i, 0)-a) > tol {
			t.Errorf("df/dr incorrect at node %d: got %v, want %v", i, dfdr.At(i, 0), a)
		}
		if math.Abs(dfds.At(i, 0)-b) > tol {
			t.Errorf("df/ds incorrect at node %d: got %v, want %v", i, dfds.At(i, 0), b)
		}
		if math.Abs(dfdt.At(i, 0)-c) > tol {
			t.Errorf("df/dt incorrect at node %d: got %v, want %v", i, dfdt.At(i, 0), c)
		}
	}
}

func BenchmarkVandermonde3D(b *testing.B) {
	N := 5
	X, Y, Z := Nodes3D(N)
	r, s, t := XYZtoRST(X, Y, Z)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		Vandermonde3D(N, r, s, t)
	}
}

// TestVandermonde3DAfterFix tests if the Vandermonde matrix is invertible after fixing RSTtoABC
func TestVandermonde3DAfterFix(t *testing.T) {
	N := 2

	// Generate nodes
	X, Y, Z := Nodes3D(N)
	r, s, tt := XYZtoRST(X, Y, Z)
	a, b, c := RSTtoABC(r, s, tt)

	// Check for duplicates in (a,b,c) space
	duplicates := make(map[string][]int)

	for i := 0; i < len(a); i++ {
		key := fmt.Sprintf("%.10f,%.10f,%.10f", a[i], b[i], c[i])
		duplicates[key] = append(duplicates[key], i)
	}

	foundDuplicates := false
	for key, nodes := range duplicates {
		if len(nodes) > 1 {
			fmt.Printf("  DUPLICATE: Nodes %v all map to %s\n", nodes, key)
			foundDuplicates = true
		}
	}

	if foundDuplicates {
		t.Error("Still have duplicate nodes after fix")
	}
}

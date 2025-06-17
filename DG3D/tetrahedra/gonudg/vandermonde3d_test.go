package gonudg

import (
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
			// Create test points
			Np := (tc.N + 1) * (tc.N + 2) * (tc.N + 3) / 6
			r, s, tt := EquidistributedNodes3D(tc.N)
			
			// Build Vandermonde matrix
			V := Vandermonde3D(tc.N, r, s, tt)
			
			// Check dimensions
			if len(V) != Np {
				t.Errorf("Wrong number of rows: got %d, want %d", len(V), Np)
			}
			if len(V[0]) != Np {
				t.Errorf("Wrong number of columns: got %d, want %d", len(V[0]), Np)
			}
			
			// Check that V is invertible (determinant != 0)
			// For small matrices, we can check condition number
			if tc.N <= 3 {
				Vinv := MatrixInverse(V)
				
				// Check V * V^{-1} = I
				I := MatrixMultiply(V, Vinv)
				for i := 0; i < Np; i++ {
					for j := 0; j < Np; j++ {
						expected := 0.0
						if i == j {
							expected = 1.0
						}
						if math.Abs(I[i][j]-expected) > 1e-10 {
							t.Errorf("V*Vinv not identity at (%d,%d): got %v", 
								i, j, I[i][j])
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
	r, s, tt := EquidistributedNodes3D(N)
	
	Vr, Vs, Vt := GradVandermonde3D(N, r, s, tt)
	
	// Check dimensions
	for _, V := range [][][]float64{Vr, Vs, Vt} {
		if len(V) != Np {
			t.Errorf("Wrong number of rows: got %d, want %d", len(V), Np)
		}
		if len(V[0]) != Np {
			t.Errorf("Wrong number of columns: got %d, want %d", len(V[0]), Np)
		}
	}
	
	// Test specific derivatives
	// For P_{1,0,0}, we expect dr/dr = const, dr/ds = function of r,s,t
	// Column for P_{1,0,0} is at index 1 (after P_{0,0,0})
	col := 1
	
	// Check that derivative values are reasonable
	for i := 0; i < Np; i++ {
		if math.IsNaN(Vr[i][col]) || math.IsNaN(Vs[i][col]) || math.IsNaN(Vt[i][col]) {
			t.Errorf("NaN in gradient at row %d", i)
		}
	}
}

func TestDmatrices3D(t *testing.T) {
	// Test differentiation matrices
	N := 2
	r, s, tt := EquidistributedNodes3D(N)
	V := Vandermonde3D(N, r, s, tt)
	
	Dr, Ds, Dt := Dmatrices3D(N, r, s, tt, V)
	
	// Test exact differentiation of polynomials
	// For a linear function f(r,s,t) = ar + bs + ct + d
	// We should get df/dr = a exactly
	
	a, b, c, d := 2.0, -3.0, 1.5, 0.5
	f := make([]float64, len(r))
	for i := 0; i < len(r); i++ {
		f[i] = a*r[i] + b*s[i] + c*tt[i] + d
	}
	
	// Compute derivatives
	dfdr := MatrixVectorMultiply(Dr, f)
	dfds := MatrixVectorMultiply(Ds, f)
	dfdt := MatrixVectorMultiply(Dt, f)
	
	// Check results
	tol := 1e-10
	for i := 0; i < len(r); i++ {
		if math.Abs(dfdr[i]-a) > tol {
			t.Errorf("df/dr incorrect at node %d: got %v, want %v", i, dfdr[i], a)
		}
		if math.Abs(dfds[i]-b) > tol {
			t.Errorf("df/ds incorrect at node %d: got %v, want %v", i, dfds[i], b)
		}
		if math.Abs(dfdt[i]-c) > tol {
			t.Errorf("df/dt incorrect at node %d: got %v, want %v", i, dfdt[i], c)
		}
	}
}

// Helper functions for testing

// EquidistributedNodes3D creates evenly spaced nodes in the reference tetrahedron
func EquidistributedNodes3D(N int) (r, s, t []float64) {
	Np := (N + 1) * (N + 2) * (N + 3) / 6
	r = make([]float64, Np)
	s = make([]float64, Np)
	t = make([]float64, Np)
	
	sk := 0
	for n := 0; n <= N; n++ {
		for m := 0; m <= N-n; m++ {
			for l := 0; l <= N-n-m; l++ {
				r[sk] = -1.0 + 2.0*float64(l)/float64(N)
				s[sk] = -1.0 + 2.0*float64(m)/float64(N)
				t[sk] = -1.0 + 2.0*float64(n)/float64(N)
				sk++
			}
		}
	}
	
	return r, s, t
}

// MatrixVectorMultiply computes y = A * x
func MatrixVectorMultiply(A [][]float64, x []float64) []float64 {
	m := len(A)
	n := len(x)
	y := make([]float64, m)
	
	for i := 0; i < m; i++ {
		sum := 0.0
		for j := 0; j < n; j++ {
			sum += A[i][j] * x[j]
		}
		y[i] = sum
	}
	
	return y
}

func BenchmarkVandermonde3D(b *testing.B) {
	N := 5
	r, s, t := EquidistributedNodes3D(N)
	
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		Vandermonde3D(N, r, s, t)
	}
}


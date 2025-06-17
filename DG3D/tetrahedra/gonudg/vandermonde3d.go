package gonudg

// Vandermonde3D initializes the 3D Vandermonde Matrix V_{ij} = phi_j(r_i, s_i, t_i)
// This is the 0-based index version of the C++ Vandermonde3D function
func Vandermonde3D(N int, r, s, t []float64) [][]float64 {
	Np := len(r)
	Ncol := (N + 1) * (N + 2) * (N + 3) / 6
	
	// Initialize the Vandermonde matrix
	V3D := make([][]float64, Np)
	for i := range V3D {
		V3D[i] = make([]float64, Ncol)
	}
	
	// Transfer to (a,b,c) coordinates
	a, b, c := RSTtoABC(r, s, t)
	
	// Build the Vandermonde matrix
	sk := 0 // 0-based column index
	for i := 0; i <= N; i++ {
		for j := 0; j <= N-i; j++ {
			for k := 0; k <= N-i-j; k++ {
				// Evaluate basis function at all points
				col := Simplex3DP(a, b, c, i, j, k)
				
				// Copy to matrix column
				for row := 0; row < Np; row++ {
					V3D[row][sk] = col[row]
				}
				sk++
			}
		}
	}
	
	return V3D
}

// GradVandermonde3D builds the gradient Vandermonde matrices
// Returns Vr, Vs, Vt where (Vr)_{ij} = dphi_j/dr at point i
func GradVandermonde3D(N int, r, s, t []float64) (Vr, Vs, Vt [][]float64) {
	Np := len(r)
	Ncol := (N + 1) * (N + 2) * (N + 3) / 6
	
	// Initialize the gradient matrices
	Vr = make([][]float64, Np)
	Vs = make([][]float64, Np)
	Vt = make([][]float64, Np)
	for i := range Vr {
		Vr[i] = make([]float64, Ncol)
		Vs[i] = make([]float64, Ncol)
		Vt[i] = make([]float64, Ncol)
	}
	
	// Build the gradient Vandermonde matrices
	sk := 0 // 0-based column index
	for i := 0; i <= N; i++ {
		for j := 0; j <= N-i; j++ {
			for k := 0; k <= N-i-j; k++ {
				// Evaluate gradient of basis function at all points
				dr, ds, dt := GradSimplex3DP(r, s, t, i, j, k)
				
				// Copy to matrix columns
				for row := 0; row < Np; row++ {
					Vr[row][sk] = dr[row]
					Vs[row][sk] = ds[row]
					Vt[row][sk] = dt[row]
				}
				sk++
			}
		}
	}
	
	return Vr, Vs, Vt
}

// Dmatrices3D computes the differentiation matrices Dr, Ds, Dt
// Given the Vandermonde matrix V and points (r,s,t)
func Dmatrices3D(N int, r, s, t []float64, V [][]float64) (Dr, Ds, Dt [][]float64) {
	// Get gradient Vandermonde matrices
	Vr, Vs, Vt := GradVandermonde3D(N, r, s, t)
	
	// Compute V inverse
	Vinv := MatrixInverse(V)
	
	// Dr = Vr * V^{-1}, etc.
	Dr = MatrixMultiply(Vr, Vinv)
	Ds = MatrixMultiply(Vs, Vinv)
	Dt = MatrixMultiply(Vt, Vinv)
	
	return Dr, Ds, Dt
}

// Helper matrix operations

// MatrixMultiply computes C = A * B
func MatrixMultiply(A, B [][]float64) [][]float64 {
	m := len(A)
	n := len(B[0])
	k := len(B)
	
	C := make([][]float64, m)
	for i := range C {
		C[i] = make([]float64, n)
	}
	
	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			sum := 0.0
			for l := 0; l < k; l++ {
				sum += A[i][l] * B[l][j]
			}
			C[i][j] = sum
		}
	}
	
	return C
}

// MatrixInverse computes the inverse of a square matrix using Gauss-Jordan elimination
func MatrixInverse(A [][]float64) [][]float64 {
	n := len(A)
	
	// Create augmented matrix [A | I]
	aug := make([][]float64, n)
	for i := range aug {
		aug[i] = make([]float64, 2*n)
		copy(aug[i], A[i])
		aug[i][n+i] = 1.0
	}
	
	// Gauss-Jordan elimination
	for i := 0; i < n; i++ {
		// Find pivot
		maxRow := i
		for k := i + 1; k < n; k++ {
			if math.Abs(aug[k][i]) > math.Abs(aug[maxRow][i]) {
				maxRow = k
			}
		}
		
		// Swap rows
		aug[i], aug[maxRow] = aug[maxRow], aug[i]
		
		// Make diagonal 1
		pivot := aug[i][i]
		for j := 0; j < 2*n; j++ {
			aug[i][j] /= pivot
		}
		
		// Eliminate column
		for k := 0; k < n; k++ {
			if k != i {
				factor := aug[k][i]
				for j := 0; j < 2*n; j++ {
					aug[k][j] -= factor * aug[i][j]
				}
			}
		}
	}
	
	// Extract inverse from augmented matrix
	inv := make([][]float64, n)
	for i := range inv {
		inv[i] = make([]float64, n)
		copy(inv[i], aug[i][n:])
	}
	
	return inv
}

// MatrixTranspose computes A^T
func MatrixTranspose(A [][]float64) [][]float64 {
	m := len(A)
	n := len(A[0])
	
	AT := make([][]float64, n)
	for i := range AT {
		AT[i] = make([]float64, m)
		for j := 0; j < m; j++ {
			AT[i][j] = A[j][i]
		}
	}
	
	return AT
}


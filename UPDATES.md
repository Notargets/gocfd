## Currently working on 3D tetrahedral FR solver

The plan is to implement tetrahedral elements based on a PKD
(Proriol-Koornwinder-Dubiner) polynomial basis. The flux will be carried on
the same order polynomials as the solution, which means Divergence and other
derivatives will require matrix multiplications across all elements. This
causes order P^6 terms in the computational cost for 3D. In order to reduce
this order to P^5, I'll be implementing the Sum Factorization approach for PKD
introduced by Fortunato and Persson in 2018.

The computational cost problem becomes clear when we look at the naive approach
to computing a derivative in 3D. Note that the below example loops are over
the DOFs of the element, which themselves number P^3, so the naive version
is order N^2 where N is O(P^3) and the Sum-factorized operation is order
N^(4/3), which substantially reduces the computational work. This becomes
critically important as we increase the order to P=4 and beyond for LES
(Large Eddy Simulation) turbulence modeling.

```
// NAIVE APPROACH: O(N^2) operations
func differentiateNaive(u []float64, dMatrix [][]float64) []float64 {
    // dMatrix is 35x35 for P=4 tetrahedron
    duDx := make([]float64, 35)
    for i := 0; i < 35; i++ {          // For each output point
        for j := 0; j < 35; j++ {      // Sum over all input points
            duDx[i] += dMatrix[i][j] * u[j]
        }
    }
    return duDx  // 35² = 1,225 operations
}
// SUM-FACTORIZED PKD: O(N^(4/3)) operations
func differentiatePKD(uModes []float64) []float64 {
    // Exploit tensor structure of PKD basis
    duDx := make([]float64, 35)
    applyDa := func(in, out []float64, rows, cols int) {
        // Apply operator along first dimension
        for j := 0; j < cols; j++ {
            for i := 0; i < rows; i++ {
                // Sum over the appropriate 1D operator
                out[i*cols+j] = /* 1D operation */
            }
        }
    }    
    // Apply three 1D operators in sequence
    temp1 := make([]float64, 35)
    applyDa(uModes, temp1, 5, 7)    // 5x7 operations
    temp2 := make([]float64, 35)
    applyDb(temp1, temp2, 5, 7)     // 5x7 operations
    applyDc(temp2, duDx)             // 35 operations
    
    return duDx  // ~105 operations total (12x faster!)
}
```

package DG1D

import "math"

func gamma0(alpha, beta float64) float64 {
    ab1 := alpha+beta + 1.
    a1 := alpha + 1.
    b1 := beta + 1.
    return math.Gamma(a1)*math.Gamma(b1)*math.Pow(2, ab1) / ab1 / math.Gamma(ab1)
}
func gamma1(alpha, beta float64) float64 {
    ab := alpha+beta
    a1 := alpha + 1.
    b1 := beta + 1.
    return a1 * b1 * gamma0(alpha, beta) / (ab+3.0)
}

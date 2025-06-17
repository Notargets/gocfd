package gonudg

import (
	"math"
)

// Simplex3DP evaluates 3D orthonormal polynomial on simplex at (a,b,c) of order (i,j,k)
// This is the 0-based index version of the C++ Simplex3DP function
func Simplex3DP(a, b, c []float64, i, j, k int) []float64 {
	n := len(a)
	P := make([]float64, n)

	// Compute Jacobi polynomials
	h1 := JacobiP(a, 0.0, 0.0, i)
	h2 := JacobiP(b, float64(2*i+1), 0.0, j)
	h3 := JacobiP(c, float64(2*(i+j)+2), 0.0, k)

	// Compute the polynomial values
	// P = 2.0*sqrt(2.0)*h1 .* h2 .* ((1.0-b).^i).*h3 .* ((1.0-c).^(i+j))
	normConst := 2.0 * math.Sqrt(2.0)

	for idx := 0; idx < n; idx++ {
		tv1 := normConst * h1[idx] * h2[idx]
		tv2 := math.Pow(1.0-b[idx], float64(i))
		tv3 := h3[idx] * math.Pow(1.0-c[idx], float64(i+j))
		P[idx] = tv1 * tv2 * tv3
	}

	return P
}

// Simplex3DPSingle evaluates 3D orthonormal polynomial at a single point
func Simplex3DPSingle(a, b, c float64, i, j, k int) float64 {
	// Compute Jacobi polynomials for single point
	h1 := JacobiPSingle(a, 0.0, 0.0, i)
	h2 := JacobiPSingle(b, float64(2*i+1), 0.0, j)
	h3 := JacobiPSingle(c, float64(2*(i+j)+2), 0.0, k)

	// Compute the polynomial value
	normConst := 2.0 * math.Sqrt(2.0)
	tv1 := normConst * h1 * h2
	tv2 := math.Pow(1.0-b, float64(i))
	tv3 := h3 * math.Pow(1.0-c, float64(i+j))

	return tv1 * tv2 * tv3
}

// GradSimplex3DP computes gradients of 3D orthonormal polynomial
// Returns the derivatives with respect to r, s, and t coordinates
func GradSimplex3DP(r, s, t []float64, id, jd, kd int) (dmodedr, dmodeds, dmodedt []float64) {
	n := len(r)
	dmodedr = make([]float64, n)
	dmodeds = make([]float64, n)
	dmodedt = make([]float64, n)

	// Convert to collapsed coordinates
	a, b, c := RSTtoABC(r, s, t)

	// Compute Jacobi polynomials and their derivatives
	fa := JacobiP(a, 0, 0, id)
	gb := JacobiP(b, float64(2*id+1), 0, jd)
	hc := JacobiP(c, float64(2*(id+jd)+2), 0, kd)

	dfa := GradJacobiP(a, 0, 0, id)
	dgb := GradJacobiP(b, float64(2*id+1), 0, jd)
	dhc := GradJacobiP(c, float64(2*(id+jd)+2), 0, kd)

	// Normalization factor
	normFactor := math.Pow(2, float64(2*id+jd)+1.5)

	// Compute each derivative component
	for i := 0; i < n; i++ {
		ai := a[i]
		bi := b[i]
		ci := c[i]

		// r-derivative
		V3Dr := dfa[i] * gb[i] * hc[i]
		if id > 0 {
			V3Dr *= math.Pow(0.5*(1.0-bi), float64(id-1))
		}
		if id+jd > 0 {
			V3Dr *= math.Pow(0.5*(1.0-ci), float64(id+jd-1))
		}

		// s-derivative
		V3Ds := 0.5 * (1.0 + ai) * V3Dr
		tmp := dgb[i] * math.Pow(0.5*(1.0-bi), float64(id))
		if id > 0 {
			tmp -= (0.5 * float64(id)) * gb[i] * math.Pow(0.5*(1.0-bi), float64(id-1))
		}
		if id+jd > 0 {
			tmp *= math.Pow(0.5*(1.0-ci), float64(id+jd-1))
		}
		tmp = fa[i] * tmp * hc[i]
		V3Ds += tmp

		// t-derivative
		V3Dt := 0.5*(1.0+ai)*V3Dr + 0.5*(1.0+bi)*tmp
		tmp2 := dhc[i] * math.Pow(0.5*(1.0-ci), float64(id+jd))
		if id+jd > 0 {
			tmp2 -= (0.5 * float64(id+jd)) * hc[i] * math.Pow(0.5*(1.0-ci), float64(id+jd-1))
		}
		tmp2 = fa[i] * gb[i] * tmp2
		tmp2 *= math.Pow(0.5*(1.0-bi), float64(id))
		V3Dt += tmp2

		// Apply normalization
		dmodedr[i] = V3Dr * normFactor
		dmodeds[i] = V3Ds * normFactor
		dmodedt[i] = V3Dt * normFactor
	}

	return
}

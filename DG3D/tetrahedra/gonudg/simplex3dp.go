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

		// Helper terms
		oneMb := 1.0 - bi
		oneMc := 1.0 - ci
		oneMbPowId := math.Pow(oneMb, float64(id))
		oneMbPowIdM1 := 1.0
		if id > 0 {
			oneMbPowIdM1 = math.Pow(oneMb, float64(id-1))
		}
		oneMcPowIJK := math.Pow(oneMc, float64(id+jd))
		oneMcPowIJKM1 := 1.0
		if id+jd > 0 {
			oneMcPowIJKM1 = math.Pow(oneMc, float64(id+jd-1))
		}

		// r-derivative
		V3Dr := dfa[i] * gb[i] * hc[i]
		if id > 0 {
			V3Dr *= oneMbPowIdM1
		}
		if id+jd > 0 {
			V3Dr *= oneMcPowIJKM1
		}
		dmodedr[i] = V3Dr * normFactor

		// s-derivative
		V3Ds := 0.5 * (1 + ai) * V3Dr
		tmp := dgb[i] * oneMbPowId
		if id > 0 {
			tmp += (-0.5 * float64(id)) * gb[i] * oneMbPowIdM1
		}
		if id+jd > 0 {
			tmp *= oneMcPowIJKM1
		}
		tmp = fa[i] * tmp * hc[i]
		V3Ds += tmp
		dmodeds[i] = V3Ds * normFactor

		// t-derivative
		V3Dt := 0.5*(1+ai)*V3Dr + 0.5*(1+bi)*tmp
		tmp2 := dhc[i] * oneMcPowIJK
		if id+jd > 0 {
			tmp2 -= 0.5 * float64(id+jd) * hc[i] * oneMcPowIJKM1
		}
		tmp2 = fa[i] * gb[i] * tmp2 * oneMbPowId
		V3Dt += tmp2
		dmodedt[i] = V3Dt * normFactor
	}

	return
}

// GradSimplex3DPSingle computes gradient at a single point
func GradSimplex3DPSingle(r, s, t float64, id, jd, kd int) (dmodedr, dmodeds, dmodedt float64) {
	// Convert to collapsed coordinates
	a, b, c := RSTtoABCSingle(r, s, t)

	// Compute Jacobi polynomials and their derivatives
	fa := JacobiPSingle(a, 0, 0, id)
	gb := JacobiPSingle(b, float64(2*id+1), 0, jd)
	hc := JacobiPSingle(c, float64(2*(id+jd)+2), 0, kd)

	dfa := GradJacobiPSingle(a, 0, 0, id)
	dgb := GradJacobiPSingle(b, float64(2*id+1), 0, jd)
	dhc := GradJacobiPSingle(c, float64(2*(id+jd)+2), 0, kd)

	// Normalization factor
	normFactor := math.Pow(2, float64(2*id+jd)+1.5)

	// Helper terms
	oneMb := 1.0 - b
	oneMc := 1.0 - c
	oneMbPowId := math.Pow(oneMb, float64(id))
	oneMbPowIdM1 := 1.0
	if id > 0 {
		oneMbPowIdM1 = math.Pow(oneMb, float64(id-1))
	}
	oneMcPowIJK := math.Pow(oneMc, float64(id+jd))
	oneMcPowIJKM1 := 1.0
	if id+jd > 0 {
		oneMcPowIJKM1 = math.Pow(oneMc, float64(id+jd-1))
	}

	// r-derivative
	V3Dr := dfa * gb * hc
	if id > 0 {
		V3Dr *= oneMbPowIdM1
	}
	if id+jd > 0 {
		V3Dr *= oneMcPowIJKM1
	}
	dmodedr = V3Dr * normFactor

	// s-derivative
	V3Ds := 0.5 * (1 + a) * V3Dr
	tmp := dgb * oneMbPowId
	if id > 0 {
		tmp += (-0.5 * float64(id)) * gb * oneMbPowIdM1
	}
	if id+jd > 0 {
		tmp *= oneMcPowIJKM1
	}
	tmp = fa * tmp * hc
	V3Ds += tmp
	dmodeds = V3Ds * normFactor

	// t-derivative
	V3Dt := 0.5*(1+a)*V3Dr + 0.5*(1+b)*tmp
	tmp2 := dhc * oneMcPowIJK
	if id+jd > 0 {
		tmp2 -= 0.5 * float64(id+jd) * hc * oneMcPowIJKM1
	}
	tmp2 = fa * gb * tmp2 * oneMbPowId
	V3Dt += tmp2
	dmodedt = V3Dt * normFactor

	return
}


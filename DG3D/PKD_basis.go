package DG3D

import (
	"github.com/notargets/gocfd/DG2D"
	"github.com/notargets/gocfd/utils"
	"math"
)

// PKDBasis represents the PKD basis for tetrahedra
type PKDBasis struct {
	P int // polynomial order
	N int // number of basis functions
}

// NewPKDBasis creates a new PKD basis of order P
func NewPKDBasis(P int) *PKDBasis {
	N := (P + 1) * (P + 2) * (P + 3) / 6
	return &PKDBasis{P: P, N: N}
}

// JacobiP computes the Jacobi polynomial of degree n with parameters alpha, beta
func JacobiP(x float64, alpha, beta float64, n int) float64 {
	if n == 0 {
		return 1.0
	}
	if n == 1 {
		return 0.5 * (alpha - beta + (alpha+beta+2)*x)
	}

	// Three-term recurrence
	p0 := 1.0
	p1 := 0.5 * (alpha - beta + (alpha+beta+2)*x)

	for k := 1; k < n; k++ {
		kF := float64(k)
		a1 := 2 * (kF + 1) * (kF + alpha + beta + 1) * (2*kF + alpha + beta)
		a2 := (2*kF + alpha + beta + 1) * (alpha*alpha - beta*beta)
		a3 := (2*kF + alpha + beta) * (2*kF + alpha + beta + 1) * (2*kF + alpha + beta + 2)
		a4 := 2 * (kF + alpha) * (kF + beta) * (2*kF + alpha + beta + 2)

		p2 := ((a2+a3*x)*p1 - a4*p0) / a1
		p0 = p1
		p1 = p2
	}

	return p1
}

// GradJacobiP computes the derivative of Jacobi polynomial
func GradJacobiP(x float64, alpha, beta float64, n int) float64 {
	if n == 0 {
		return 0.0
	}
	return 0.5 * (alpha + beta + float64(n) + 1) * JacobiP(x, alpha+1, beta+1, n-1)
}

// EvalBasis evaluates the PKD basis at a single point (r,s,t)
func (basis *PKDBasis) EvalBasis(r, s, t float64) []float64 {
	phi := make([]float64, basis.N)

	// Transform to orthogonal coordinates
	eps := 1e-10

	// Handle boundary cases first
	var a, b, c float64
	if math.Abs(1-s-t) < eps {
		a = -1
	} else {
		a = 2*(1+r)/(1-s-t) - 1
	}

	if math.Abs(1-t) < eps {
		b = -1
	} else {
		b = 2*(1+s)/(1-t) - 1
	}

	c = 2*(1+t) - 1

	idx := 0
	for i := 0; i <= basis.P; i++ {
		for j := 0; j <= basis.P-i; j++ {
			for k := 0; k <= basis.P-i-j; k++ {
				// Evaluate PKD basis function
				val := JacobiP(a, 0, 0, i)
				val *= JacobiP(b, float64(2*i+1), 0, j)
				val *= JacobiP(c, float64(2*i+2*j+2), 0, k)

				// Scaling factors - handle edge cases
				if i > 0 {
					val *= math.Pow((1-b)/2, float64(i))
				}
				if i+j > 0 {
					val *= math.Pow((1-c)/2, float64(i+j))
				}

				// Overall scaling
				val *= math.Pow(2, float64(i+j+k))

				phi[idx] = val
				idx++
			}
		}
	}

	return phi
}

// Vandermonde generates the Vandermonde matrix for given evaluation points
func (basis *PKDBasis) Vandermonde(r, s, t []float64) utils.Matrix {
	Np := len(r)
	V := utils.NewMatrix(Np, basis.N)

	for i := 0; i < Np; i++ {
		phi := basis.EvalBasis(r[i], s[i], t[i])
		for j := 0; j < basis.N; j++ {
			V.Set(i, j, phi[j])
		}
	}

	return V
}

// DerivativeMatrix computes the derivative matrix in direction dir (0=r, 1=s, 2=t)
func (basis *PKDBasis) DerivativeMatrix(r, s, t []float64, dir int) utils.Matrix {
	Np := len(r)
	Dr := utils.NewMatrix(Np, basis.N)

	for i := 0; i < Np; i++ {
		dphi := basis.EvalDerivBasis(r[i], s[i], t[i], dir)
		for j := 0; j < basis.N; j++ {
			Dr.Set(i, j, dphi[j])
		}
	}

	// Compute Dr = Dr * inv(V)
	V := basis.Vandermonde(r, s, t)

	// Check condition number
	cond := V.ConditionNumber()
	if cond > 1e12 {
		// If badly conditioned, try pseudoinverse
		Vinv := V.PseudoInverse(1e-10)
		return Dr.Mul(Vinv)
	}

	Vinv, err := V.Inverse()
	if err != nil {
		// Try pseudoinverse as fallback
		Vinv = V.PseudoInverse(1e-10)
	}

	return Dr.Mul(Vinv)
}

// EvalDerivBasis evaluates derivatives of PKD basis at a single point
func (basis *PKDBasis) EvalDerivBasis(r, s, t float64, dir int) []float64 {
	dphi := make([]float64, basis.N)

	// Transform to orthogonal coordinates with small epsilon to avoid division by zero
	eps := 1e-10
	a := 2*(1+r)/(1-s-t+eps) - 1
	b := 2*(1+s)/(1-t+eps) - 1
	c := 2*(1+t) - 1

	// Handle boundary cases
	if math.Abs(s+t-1) < eps {
		a = -1
	}
	if math.Abs(t-1) < eps {
		b = -1
		c = 1
	}

	// Compute coordinate derivatives
	var da_dr, da_ds, da_dt float64
	var db_dr, db_ds, db_dt float64
	var dc_dr, dc_ds, dc_dt float64

	if math.Abs(s+t-1) > eps {
		denom := 1 - s - t + eps
		da_dr = 2 / denom
		da_ds = 2 * (1 + r) / (denom * denom)
		da_dt = 2 * (1 + r) / (denom * denom)
	}

	if math.Abs(t-1) > eps {
		denom := 1 - t + eps
		db_dr = 0
		db_ds = 2 / denom
		db_dt = 2 * (1 + s) / (denom * denom)
	}

	dc_dr = 0
	dc_ds = 0
	dc_dt = 2

	idx := 0
	for i := 0; i <= basis.P; i++ {
		for j := 0; j <= basis.P-i; j++ {
			for k := 0; k <= basis.P-i-j; k++ {
				// Compute basis function value and derivatives
				Pa := JacobiP(a, 0, 0, i)
				Pb := JacobiP(b, float64(2*i+1), 0, j)
				Pc := JacobiP(c, float64(2*i+2*j+2), 0, k)

				dPa := GradJacobiP(a, 0, 0, i)
				dPb := GradJacobiP(b, float64(2*i+1), 0, j)
				dPc := GradJacobiP(c, float64(2*i+2*j+2), 0, k)

				// Scaling factors (with protection against 0^0)
				scale1 := 1.0
				if i > 0 || (1-b)/2 != 0 {
					scale1 = math.Pow((1-b)/2, float64(i))
				}
				scale2 := 1.0
				if i+j > 0 || (1-c)/2 != 0 {
					scale2 = math.Pow((1-c)/2, float64(i+j))
				}
				scale3 := math.Pow(2, float64(i+j+k))

				// Chain rule for derivatives
				var dphi_dr, dphi_ds, dphi_dt float64

				// Derivative with respect to a
				term1 := dPa * Pb * Pc * scale1 * scale2 * scale3
				dphi_dr += term1 * da_dr
				dphi_ds += term1 * da_ds
				dphi_dt += term1 * da_dt

				// Derivative with respect to b (including scale derivative)
				term2 := Pa * dPb * Pc * scale1 * scale2 * scale3
				if i > 0 && math.Abs(1-b) > eps {
					dscale1_db := float64(i) * math.Pow((1-b)/2, float64(i-1)) * (-0.5)
					term2 += Pa * Pb * Pc * dscale1_db * scale2 * scale3
				}
				dphi_dr += term2 * db_dr
				dphi_ds += term2 * db_ds
				dphi_dt += term2 * db_dt

				// Derivative with respect to c (including scale derivative)
				term3 := Pa * Pb * dPc * scale1 * scale2 * scale3
				if i+j > 0 && math.Abs(1-c) > eps {
					dscale2_dc := float64(i+j) * math.Pow((1-c)/2, float64(i+j-1)) * (-0.5)
					term3 += Pa * Pb * Pc * scale1 * dscale2_dc * scale3
				}
				dphi_dr += term3 * dc_dr
				dphi_ds += term3 * dc_ds
				dphi_dt += term3 * dc_dt

				// Store the appropriate derivative
				switch dir {
				case 0:
					dphi[idx] = dphi_dr
				case 1:
					dphi[idx] = dphi_ds
				case 2:
					dphi[idx] = dphi_dt
				}
				idx++
			}
		}
	}

	return dphi
}

// SumFactorization implements the sum factorization operators
// Based on the collapsed coordinate approach for tetrahedra
type SumFactorization struct {
	P   int
	N   int
	n1d int // 1D basis size
	n2d int // 2D triangle basis size

	// 1D operators for the tensor product structure
	V1d utils.Matrix // 1D Vandermonde
	D1d utils.Matrix // 1D derivative

	// 2D operators for triangular faces
	V2d  utils.Matrix // 2D Vandermonde on triangle
	Dr2d utils.Matrix // 2D r-derivative on triangle
	Ds2d utils.Matrix // 2D s-derivative on triangle

	// Quadrature weights for collapsed coordinates
	wa []float64 // weights for 'a' coordinate
	wb []float64 // weights for 'b' coordinate
	wc []float64 // weights for 'c' coordinate

	// Store nodes for quadrature
	nodes1d []float64
	nodes2d struct {
		r []float64
		s []float64
	}
}

// NewSumFactorization creates sum factorized operators using collapsed coordinates
func NewSumFactorization(P int) *SumFactorization {
	N := (P + 1) * (P + 2) * (P + 3) / 6
	n1d := P + 1
	n2d := (P + 1) * (P + 2) / 2

	sf := &SumFactorization{
		P:   P,
		N:   N,
		n1d: n1d,
		n2d: n2d,
	}

	// Initialize 1D operators (Legendre basis on [-1,1])
	sf.init1DOperators()

	// Initialize 2D operators (PKD basis on triangle)
	sf.init2DOperators()

	// Initialize quadrature weights for collapsed coordinates
	sf.initQuadratureWeights()

	return sf
}

// init1DOperators initializes 1D Legendre operators
func (sf *SumFactorization) init1DOperators() {
	// 1D Legendre-Gauss-Lobatto nodes
	sf.nodes1d = make([]float64, sf.n1d)
	for i := 0; i < sf.n1d; i++ {
		sf.nodes1d[i] = -math.Cos(math.Pi * float64(i) / float64(sf.P))
	}

	// 1D Vandermonde matrix
	sf.V1d = utils.NewMatrix(sf.n1d, sf.n1d)
	for i := 0; i < sf.n1d; i++ {
		for j := 0; j < sf.n1d; j++ {
			sf.V1d.Set(i, j, JacobiP(sf.nodes1d[i], 0, 0, j))
		}
	}

	// 1D derivative matrix
	sf.D1d = utils.NewMatrix(sf.n1d, sf.n1d)
	for i := 0; i < sf.n1d; i++ {
		for j := 0; j < sf.n1d; j++ {
			sf.D1d.Set(i, j, GradJacobiP(sf.nodes1d[i], 0, 0, j))
		}
	}

	// D1d = D1d * inv(V1d)
	Vinv, err := sf.V1d.Inverse()
	if err != nil {
		panic(err)
	}
	sf.D1d = sf.D1d.Mul(Vinv)
}

// init2DOperators initializes 2D operators on reference triangle
func (sf *SumFactorization) init2DOperators() {
	// Generate 2D nodes on reference triangle using equispaced distribution
	sf.nodes2d.r = make([]float64, sf.n2d)
	sf.nodes2d.s = make([]float64, sf.n2d)

	idx := 0
	for i := 0; i <= sf.P; i++ {
		for j := 0; j <= sf.P-i; j++ {
			// Map to reference triangle [-1,1] x [-1,1-r]
			sf.nodes2d.r[idx] = -1.0 + 2.0*float64(j)/float64(sf.P)
			sf.nodes2d.s[idx] = -1.0 + 2.0*float64(i)/float64(sf.P)
			idx++
		}
	}

	// Create 2D basis
	r2d := utils.NewVector(sf.n2d, sf.nodes2d.r)
	s2d := utils.NewVector(sf.n2d, sf.nodes2d.s)

	// Use DG2D Jacobi basis
	jb2d := DG2D.NewJacobiBasis2D(sf.P, r2d, s2d, 0.0, 0.0)

	// Get Vandermonde and derivative matrices
	sf.V2d = jb2d.V
	sf.Dr2d = jb2d.Vr.Mul(jb2d.Vinv)
	sf.Ds2d = jb2d.Vs.Mul(jb2d.Vinv)
}

// initQuadratureWeights initializes weights for collapsed coordinate quadrature
func (sf *SumFactorization) initQuadratureWeights() {
	// Gauss-Jacobi weights for the collapsed coordinates
	// These account for the Jacobian of the transformation
	sf.wa = make([]float64, sf.n1d)
	sf.wb = make([]float64, sf.n1d)
	sf.wc = make([]float64, sf.n1d)

	// For collapsed coordinates, we need Gauss-Jacobi quadrature with specific parameters
	// a: Gauss-Legendre (alpha=0, beta=0)
	// b: Gauss-Jacobi (alpha=1, beta=0)
	// c: Gauss-Jacobi (alpha=2, beta=0)

	// Get Gauss-Legendre weights for 'a' coordinate
	for i := 0; i < sf.n1d; i++ {
		// Gauss-Lobatto weights approximation
		if i == 0 || i == sf.n1d-1 {
			sf.wa[i] = 2.0 / float64(sf.n1d*(sf.n1d-1))
		} else {
			// Interior points - use Legendre polynomial derivative
			x := sf.nodes1d[i]
			Pn := JacobiP(x, 0, 0, sf.P)
			sf.wa[i] = 2.0 / float64(sf.n1d*(sf.n1d-1)) / (Pn * Pn)
		}
	}

	// Weights for 'b' coordinate include (1-b)^1 factor
	for i := 0; i < sf.n1d; i++ {
		sf.wb[i] = sf.wa[i] * math.Pow(1.0+sf.nodes1d[i], 1.0) / 2.0
	}

	// Weights for 'c' coordinate include (1-c)^2 factor
	for i := 0; i < sf.n1d; i++ {
		sf.wc[i] = sf.wa[i] * math.Pow(1.0+sf.nodes1d[i], 2.0) / 4.0
	}
}

// ReshapeToTensor reshapes tetrahedral data to tensor format for sum factorization
func (sf *SumFactorization) ReshapeToTensor(u []float64) [][][]float64 {
	// Map from tetrahedral PKD basis ordering to tensor product basis
	// using the collapsed coordinate transformation
	tensor := make([][][]float64, sf.n1d)
	for i := range tensor {
		tensor[i] = make([][]float64, sf.n1d)
		for j := range tensor[i] {
			tensor[i][j] = make([]float64, sf.n1d)
		}
	}

	// PKD basis ordering: loop over i+j+k <= P
	// We need to map this to tensor indices (i,j,k) where each runs 0 to P
	// In collapsed coordinates, we still have the constraint that modes are zero
	// when indices exceed the polynomial order

	idx := 0
	for i := 0; i <= sf.P; i++ {
		for j := 0; j <= sf.P-i; j++ {
			for k := 0; k <= sf.P-i-j; k++ {
				if idx < len(u) && i < sf.n1d && j < sf.n1d && k < sf.n1d {
					tensor[i][j][k] = u[idx]
				}
				idx++
			}
		}
	}

	// Zero out tensor entries that don't correspond to tetrahedral modes
	for i := 0; i < sf.n1d; i++ {
		for j := 0; j < sf.n1d; j++ {
			for k := 0; k < sf.n1d; k++ {
				if i+j+k > sf.P {
					tensor[i][j][k] = 0.0
				}
			}
		}
	}

	return tensor
}

// ReshapeFromTensor converts tensor format back to tetrahedral basis
func (sf *SumFactorization) ReshapeFromTensor(tensor [][][]float64) []float64 {
	u := make([]float64, sf.N)

	idx := 0
	for i := 0; i <= sf.P; i++ {
		for j := 0; j <= sf.P-i; j++ {
			for k := 0; k <= sf.P-i-j; k++ {
				if i < sf.n1d && j < sf.n1d && k < sf.n1d {
					u[idx] = tensor[i][j][k]
				}
				idx++
			}
		}
	}

	return u
}

// ApplyDr applies the r-derivative using sum factorization
func (sf *SumFactorization) ApplyDr(u []float64) []float64 {
	// Convert to tensor format
	U := sf.ReshapeToTensor(u)
	DU := make([][][]float64, sf.n1d)
	for i := range DU {
		DU[i] = make([][]float64, sf.n1d)
		for j := range DU[i] {
			DU[i][j] = make([]float64, sf.n1d)
		}
	}

	// Apply tensor product of operators
	// In collapsed coordinates: d/dr affects all three coordinates
	// d/dr = (2/(1-s-t)) d/da

	// First apply d/da
	for j := 0; j < sf.n1d; j++ {
		for k := 0; k < sf.n1d; k++ {
			// Only process non-zero modes
			if j+k <= sf.P {
				// Extract 1D slice
				slice := make([]float64, sf.n1d)
				for i := 0; i < sf.n1d && i <= sf.P-j-k; i++ {
					slice[i] = U[i][j][k]
				}

				// Apply 1D derivative
				dslice := sf.apply1DDeriv(slice)

				// Store result
				for i := 0; i < sf.n1d && i <= sf.P-j-k; i++ {
					DU[i][j][k] = dslice[i]
				}
			}
		}
	}

	// Transform derivatives to physical coordinates
	sf.transformDerivativesFromCollapsed(DU, 0)

	// Convert back to vector format
	return sf.ReshapeFromTensor(DU)
}

// ApplyDs applies the s-derivative using sum factorization
func (sf *SumFactorization) ApplyDs(u []float64) []float64 {
	U := sf.ReshapeToTensor(u)
	DU := make([][][]float64, sf.n1d)
	for i := range DU {
		DU[i] = make([][]float64, sf.n1d)
		for j := range DU[i] {
			DU[i][j] = make([]float64, sf.n1d)
		}
	}

	// In collapsed coordinates: d/ds = ((1+a)/(1-t)) d/da + (2/(1-t)) d/db
	// This requires two tensor operations

	// First term: ((1+a)/(1-t)) d/da
	DU1 := sf.applyDaForDs(U)

	// Second term: (2/(1-t)) d/db
	DU2 := sf.applyDbForDs(U)

	// Combine
	for i := 0; i < sf.n1d; i++ {
		for j := 0; j < sf.n1d; j++ {
			for k := 0; k < sf.n1d; k++ {
				DU[i][j][k] = DU1[i][j][k] + DU2[i][j][k]
			}
		}
	}

	sf.transformDerivativesFromCollapsed(DU, 1)
	return sf.ReshapeFromTensor(DU)
}

// ApplyDt applies the t-derivative using sum factorization
func (sf *SumFactorization) ApplyDt(u []float64) []float64 {
	U := sf.ReshapeToTensor(u)
	DU := make([][][]float64, sf.n1d)
	for i := range DU {
		DU[i] = make([][]float64, sf.n1d)
		for j := range DU[i] {
			DU[i][j] = make([]float64, sf.n1d)
		}
	}

	// In collapsed coordinates: d/dt involves all three coordinates
	// d/dt = ((1+a)/(1-s-t)) d/da + ((1+b)/(1-t)) d/db + 2 d/dc

	// Apply the three terms
	DU1 := sf.applyDaForDt(U)
	DU2 := sf.applyDbForDt(U)
	DU3 := sf.applyDc(U)

	// Combine
	for i := 0; i < sf.n1d; i++ {
		for j := 0; j < sf.n1d; j++ {
			for k := 0; k < sf.n1d; k++ {
				DU[i][j][k] = DU1[i][j][k] + DU2[i][j][k] + DU3[i][j][k]
			}
		}
	}

	sf.transformDerivativesFromCollapsed(DU, 2)
	return sf.ReshapeFromTensor(DU)
}

// Helper functions for applying derivatives in collapsed coordinates
func (sf *SumFactorization) applyDaForDs(U [][][]float64) [][][]float64 {
	DU := make([][][]float64, sf.n1d)
	for i := range DU {
		DU[i] = make([][]float64, sf.n1d)
		for j := range DU[i] {
			DU[i][j] = make([]float64, sf.n1d)
		}
	}

	// Apply d/da with scaling (1+a)
	for j := 0; j < sf.n1d; j++ {
		for k := 0; k < sf.n1d; k++ {
			if j+k <= sf.P {
				slice := make([]float64, sf.n1d)
				for i := 0; i < sf.n1d && i <= sf.P-j-k; i++ {
					slice[i] = U[i][j][k]
				}

				dslice := sf.apply1DDeriv(slice)

				for i := 0; i < sf.n1d && i <= sf.P-j-k; i++ {
					// Scale by (1+a) where a = nodes1d[i]
					DU[i][j][k] = dslice[i] * (1.0 + sf.nodes1d[i])
				}
			}
		}
	}

	return DU
}

func (sf *SumFactorization) applyDbForDs(U [][][]float64) [][][]float64 {
	DU := make([][][]float64, sf.n1d)
	for i := range DU {
		DU[i] = make([][]float64, sf.n1d)
		for j := range DU[i] {
			DU[i][j] = make([]float64, sf.n1d)
		}
	}

	// Apply d/db (second index)
	for i := 0; i < sf.n1d; i++ {
		for k := 0; k < sf.n1d; k++ {
			if i+k <= sf.P {
				slice := make([]float64, sf.n1d)
				for j := 0; j < sf.n1d && i+j+k <= sf.P; j++ {
					slice[j] = U[i][j][k]
				}

				dslice := sf.apply1DDeriv(slice)

				for j := 0; j < sf.n1d && i+j+k <= sf.P; j++ {
					DU[i][j][k] = 2.0 * dslice[j]
				}
			}
		}
	}

	return DU
}

func (sf *SumFactorization) applyDaForDt(U [][][]float64) [][][]float64 {
	DU := make([][][]float64, sf.n1d)
	for i := range DU {
		DU[i] = make([][]float64, sf.n1d)
		for j := range DU[i] {
			DU[i][j] = make([]float64, sf.n1d)
		}
	}

	// Apply d/da with scaling (1+a)
	for j := 0; j < sf.n1d; j++ {
		for k := 0; k < sf.n1d; k++ {
			if j+k <= sf.P {
				slice := make([]float64, sf.n1d)
				for i := 0; i < sf.n1d && i <= sf.P-j-k; i++ {
					slice[i] = U[i][j][k]
				}

				dslice := sf.apply1DDeriv(slice)

				for i := 0; i < sf.n1d && i <= sf.P-j-k; i++ {
					DU[i][j][k] = dslice[i] * (1.0 + sf.nodes1d[i])
				}
			}
		}
	}

	return DU
}

func (sf *SumFactorization) applyDbForDt(U [][][]float64) [][][]float64 {
	DU := make([][][]float64, sf.n1d)
	for i := range DU {
		DU[i] = make([][]float64, sf.n1d)
		for j := range DU[i] {
			DU[i][j] = make([]float64, sf.n1d)
		}
	}

	// Apply d/db with scaling (1+b)
	for i := 0; i < sf.n1d; i++ {
		for k := 0; k < sf.n1d; k++ {
			if i+k <= sf.P {
				slice := make([]float64, sf.n1d)
				for j := 0; j < sf.n1d && i+j+k <= sf.P; j++ {
					slice[j] = U[i][j][k]
				}

				dslice := sf.apply1DDeriv(slice)

				for j := 0; j < sf.n1d && i+j+k <= sf.P; j++ {
					DU[i][j][k] = dslice[j] * (1.0 + sf.nodes1d[j])
				}
			}
		}
	}

	return DU
}

func (sf *SumFactorization) applyDc(U [][][]float64) [][][]float64 {
	DU := make([][][]float64, sf.n1d)
	for i := range DU {
		DU[i] = make([][]float64, sf.n1d)
		for j := range DU[i] {
			DU[i][j] = make([]float64, sf.n1d)
		}
	}

	// Apply d/dc (third index)
	for i := 0; i < sf.n1d; i++ {
		for j := 0; j < sf.n1d; j++ {
			if i+j <= sf.P {
				slice := make([]float64, sf.n1d)
				for k := 0; k < sf.n1d && i+j+k <= sf.P; k++ {
					slice[k] = U[i][j][k]
				}

				dslice := sf.apply1DDeriv(slice)

				for k := 0; k < sf.n1d && i+j+k <= sf.P; k++ {
					DU[i][j][k] = 2.0 * dslice[k]
				}
			}
		}
	}

	return DU
}

// apply1DDeriv applies the 1D derivative operator to a vector
func (sf *SumFactorization) apply1DDeriv(u []float64) []float64 {
	du := make([]float64, len(u))
	for i := 0; i < len(u); i++ {
		for j := 0; j < len(u); j++ {
			du[i] += sf.D1d.At(i, j) * u[j]
		}
	}
	return du
}

// transformDerivativesFromCollapsed accounts for the chain rule from collapsed coordinates
func (sf *SumFactorization) transformDerivativesFromCollapsed(DU [][][]float64, dir int) {
	// This implements the chain rule for the collapsed coordinate transformation
	// (a,b,c) -> (r,s,t) where:
	// a = 2(1+r)/(1-s-t) - 1
	// b = 2(1+s)/(1-t) - 1
	// c = 2(1+t) - 1
	//
	// The derivatives have already been computed in collapsed coordinates.
	// Now we need to apply the proper scaling factors based on the chain rule.

	// The scaling has already been partially applied in the Apply* methods
	// This function applies any remaining coordinate-dependent factors

	for i := 0; i < sf.n1d; i++ {
		for j := 0; j < sf.n1d; j++ {
			for k := 0; k < sf.n1d; k++ {
				if i+j+k <= sf.P {
					// Convert indices to coordinates
					a := sf.nodes1d[i]
					b := sf.nodes1d[j]
					c := sf.nodes1d[k]
					_ = a

					// Compute scaling factors based on collapsed coordinate Jacobian
					switch dir {
					case 0: // r-derivative
						// d/dr has factor 2/(1-s-t) = 2/((1-b)(1-c)/2)
						scale := 4.0 / ((1.0-b)*(1.0-c) + 1e-10)
						DU[i][j][k] *= scale
					case 1: // s-derivative
						// d/ds has factor 2/(1-t) = 2/((1-c)/2)
						scale := 4.0 / (1.0 - c + 1e-10)
						DU[i][j][k] *= scale
					case 2: // t-derivative
						// d/dt already has correct scaling
					}
				}
			}
		}
	}
}

// AffineTransform handles the transformation from reference to physical coordinates
type AffineTransform struct {
	// Jacobian components
	xr, xs, xt float64
	yr, ys, yt float64
	zr, zs, zt float64
	J          float64 // Jacobian determinant
}

// NewAffineTransform creates an affine transformation from tetrahedral vertices
func NewAffineTransform(v0, v1, v2, v3 [3]float64) *AffineTransform {
	at := &AffineTransform{}

	// Compute Jacobian matrix
	at.xr = v1[0] - v0[0]
	at.xs = v2[0] - v0[0]
	at.xt = v3[0] - v0[0]

	at.yr = v1[1] - v0[1]
	at.ys = v2[1] - v0[1]
	at.yt = v3[1] - v0[1]

	at.zr = v1[2] - v0[2]
	at.zs = v2[2] - v0[2]
	at.zt = v3[2] - v0[2]

	// Compute Jacobian determinant
	at.J = at.xr*(at.ys*at.zt-at.yt*at.zs) -
		at.xs*(at.yr*at.zt-at.yt*at.zr) +
		at.xt*(at.yr*at.zs-at.ys*at.zr)

	return at
}

// TransformDerivatives transforms derivatives from reference to physical coordinates
func (at *AffineTransform) TransformDerivatives(dr, ds, dt []float64) (dx, dy, dz []float64) {
	n := len(dr)
	dx = make([]float64, n)
	dy = make([]float64, n)
	dz = make([]float64, n)

	// Compute inverse Jacobian components
	invJ := 1.0 / at.J

	// Inverse Jacobian matrix components
	rx := invJ * (at.ys*at.zt - at.yt*at.zs)
	ry := invJ * (at.yt*at.xs - at.ys*at.xt)
	rz := invJ * (at.xs*at.zt - at.xt*at.zs)

	sx := invJ * (at.zr*at.yt - at.yr*at.zt)
	sy := invJ * (at.xr*at.zt - at.xt*at.zr)
	sz := invJ * (at.yr*at.xt - at.xr*at.yt)

	tx := invJ * (at.yr*at.zs - at.zr*at.ys)
	ty := invJ * (at.zr*at.xs - at.xr*at.zs)
	tz := invJ * (at.xr*at.ys - at.yr*at.xs)

	// Apply transformation
	for i := 0; i < n; i++ {
		dx[i] = rx*dr[i] + sx*ds[i] + tx*dt[i]
		dy[i] = ry*dr[i] + sy*ds[i] + ty*dt[i]
		dz[i] = rz*dr[i] + sz*ds[i] + tz*dt[i]
	}

	return dx, dy, dz
}

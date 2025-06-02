package DG3D

import (
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
	a := 2*(1+r)/(1-s-t) - 1
	b := 2*(1+s)/(1-t) - 1
	c := 2*(1+t) - 1

	if math.Abs(s+t-1) < 1e-15 {
		a = -1
	}
	if math.Abs(t-1) < 1e-15 {
		b = -1
	}

	idx := 0
	for i := 0; i <= basis.P; i++ {
		for j := 0; j <= basis.P-i; j++ {
			for k := 0; k <= basis.P-i-j; k++ {
				// Evaluate PKD basis function
				val := JacobiP(a, 0, 0, i)
				val *= JacobiP(b, float64(2*i+1), 0, j)
				val *= JacobiP(c, float64(2*i+2*j+2), 0, k)

				// Scaling factors
				val *= math.Pow((1-b)/2, float64(i))
				val *= math.Pow((1-c)/2, float64(i+j))
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
	Vinv, err := V.Inverse()
	if err != nil {
		panic(err)
	}

	return Dr.Mul(Vinv)
}

// EvalDerivBasis evaluates derivatives of PKD basis at a single point
func (basis *PKDBasis) EvalDerivBasis(r, s, t float64, dir int) []float64 {
	dphi := make([]float64, basis.N)

	// Transform to orthogonal coordinates
	a := 2*(1+r)/(1-s-t) - 1
	b := 2*(1+s)/(1-t) - 1
	c := 2*(1+t) - 1

	if math.Abs(s+t-1) < 1e-15 {
		a = -1
	}
	if math.Abs(t-1) < 1e-15 {
		b = -1
	}

	// Compute coordinate derivatives
	var da_dr, da_ds, da_dt float64
	var db_dr, db_ds, db_dt float64
	var dc_dr, dc_ds, dc_dt float64

	if math.Abs(s+t-1) > 1e-15 {
		da_dr = 2 / (1 - s - t)
		da_ds = 2 * (1 + r) / ((1 - s - t) * (1 - s - t))
		da_dt = 2 * (1 + r) / ((1 - s - t) * (1 - s - t))
	}

	if math.Abs(t-1) > 1e-15 {
		db_dr = 0
		db_ds = 2 / (1 - t)
		db_dt = 2 * (1 + s) / ((1 - t) * (1 - t))
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

				// Scaling factors
				scale1 := math.Pow((1-b)/2, float64(i))
				scale2 := math.Pow((1-c)/2, float64(i+j))
				scale3 := math.Pow(2, float64(i+j+k))

				// Chain rule for derivatives
				var dphi_dr, dphi_ds, dphi_dt float64

				// Derivative with respect to a
				term1 := dPa * Pb * Pc * scale1 * scale2 * scale3
				dphi_dr += term1 * da_dr
				dphi_ds += term1 * da_ds
				dphi_dt += term1 * da_dt

				// Derivative with respect to b
				term2 := Pa * dPb * Pc * scale1 * scale2 * scale3
				if i > 0 {
					term2 += Pa * Pb * Pc * float64(i) * math.Pow((1-b)/2, float64(i-1)) * (-0.5) * scale2 * scale3
				}
				dphi_dr += term2 * db_dr
				dphi_ds += term2 * db_ds
				dphi_dt += term2 * db_dt

				// Derivative with respect to c
				term3 := Pa * Pb * dPc * scale1 * scale2 * scale3
				if i+j > 0 {
					term3 += Pa * Pb * Pc * scale1 * float64(i+j) * math.Pow((1-c)/2, float64(i+j-1)) * (-0.5) * scale3
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
	wa []float64
	wb []float64
	wc []float64
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
	nodes1d := make([]float64, sf.n1d)
	for i := 0; i < sf.n1d; i++ {
		nodes1d[i] = -math.Cos(math.Pi * float64(i) / float64(sf.P))
	}

	// 1D Vandermonde matrix
	sf.V1d = utils.NewMatrix(sf.n1d, sf.n1d)
	for i := 0; i < sf.n1d; i++ {
		for j := 0; j < sf.n1d; j++ {
			sf.V1d.Set(i, j, JacobiP(nodes1d[i], 0, 0, j))
		}
	}

	// 1D derivative matrix
	sf.D1d = utils.NewMatrix(sf.n1d, sf.n1d)
	for i := 0; i < sf.n1d; i++ {
		for j := 0; j < sf.n1d; j++ {
			sf.D1d.Set(i, j, GradJacobiP(nodes1d[i], 0, 0, j))
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
	// Generate 2D nodes on reference triangle
	r2d := make([]float64, sf.n2d)
	s2d := make([]float64, sf.n2d)

	idx := 0
	for i := 0; i <= sf.P; i++ {
		for j := 0; j <= sf.P-i; j++ {
			r2d[idx] = -1.0 + 2.0*float64(i)/float64(sf.P)
			s2d[idx] = -1.0 + 2.0*float64(j)/float64(sf.P)
			idx++
		}
	}

	// Build 2D PKD basis on triangle (for P=P, N=n2d, no t-coordinate)
	// We need a 2D version - simplified here
	sf.V2d = utils.NewMatrix(sf.n2d, sf.n2d)
	sf.Dr2d = utils.NewMatrix(sf.n2d, sf.n2d)
	sf.Ds2d = utils.NewMatrix(sf.n2d, sf.n2d)

	// Initialize with identity for now - full implementation would use 2D PKD basis
	for i := 0; i < sf.n2d; i++ {
		sf.V2d.Set(i, i, 1.0)
		sf.Dr2d.Set(i, i, 0.0)
		sf.Ds2d.Set(i, i, 0.0)
	}
}

// initQuadratureWeights initializes weights for collapsed coordinate quadrature
func (sf *SumFactorization) initQuadratureWeights() {
	// Weights for collapsed coordinate transformation
	// These account for the Jacobian of the transformation
	sf.wa = make([]float64, sf.n1d)
	sf.wb = make([]float64, sf.n1d)
	sf.wc = make([]float64, sf.n1d)

	// Gauss-Jacobi weights for the collapsed coordinates
	for i := 0; i < sf.n1d; i++ {
		x := -math.Cos(math.Pi * float64(i) / float64(sf.P))
		// Weight includes Jacobian factors from collapsed transformation
		sf.wa[i] = 2.0 / float64(sf.n1d) // Simplified - actual weights more complex
		sf.wb[i] = math.Pow(1.0-x, 1.0) * sf.wa[i]
		sf.wc[i] = math.Pow(1.0-x, 2.0) * sf.wa[i]
	}
}

// ReshapeToTensor reshapes tetrahedral data to tensor format for sum factorization
func (sf *SumFactorization) ReshapeToTensor(u []float64) [][][]float64 {
	// Map from tetrahedral basis to tensor product basis
	// using the collapsed coordinate transformation
	tensor := make([][][]float64, sf.n1d)
	for i := range tensor {
		tensor[i] = make([][]float64, sf.n1d)
		for j := range tensor[i] {
			tensor[i][j] = make([]float64, sf.n1d)
		}
	}

	// This is a simplified mapping - actual implementation requires
	// careful indexing based on the PKD basis ordering
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
	// Dr = D1d ⊗ I ⊗ I (in collapsed coordinates)
	for j := 0; j < sf.n1d; j++ {
		for k := 0; k < sf.n1d; k++ {
			// Extract 1D slice
			slice := make([]float64, sf.n1d)
			for i := 0; i < sf.n1d; i++ {
				slice[i] = U[i][j][k]
			}

			// Apply 1D derivative
			dslice := sf.apply1DDeriv(slice)

			// Store result
			for i := 0; i < sf.n1d; i++ {
				DU[i][j][k] = dslice[i]
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

	// Apply tensor product: Ds = I ⊗ D1d ⊗ I
	for i := 0; i < sf.n1d; i++ {
		for k := 0; k < sf.n1d; k++ {
			slice := make([]float64, sf.n1d)
			for j := 0; j < sf.n1d; j++ {
				slice[j] = U[i][j][k]
			}

			dslice := sf.apply1DDeriv(slice)

			for j := 0; j < sf.n1d; j++ {
				DU[i][j][k] = dslice[j]
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

	// Apply tensor product: Dt = I ⊗ I ⊗ D1d
	for i := 0; i < sf.n1d; i++ {
		for j := 0; j < sf.n1d; j++ {
			slice := make([]float64, sf.n1d)
			for k := 0; k < sf.n1d; k++ {
				slice[k] = U[i][j][k]
			}

			dslice := sf.apply1DDeriv(slice)

			for k := 0; k < sf.n1d; k++ {
				DU[i][j][k] = dslice[k]
			}
		}
	}

	sf.transformDerivativesFromCollapsed(DU, 2)
	return sf.ReshapeFromTensor(DU)
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

	// The actual implementation would apply the appropriate Jacobian factors
	// This is simplified for demonstration
	for i := 0; i < sf.n1d; i++ {
		for j := 0; j < sf.n1d; j++ {
			for k := 0; k < sf.n1d; k++ {
				// Apply coordinate transformation scaling
				switch dir {
				case 0: // r-derivative
					DU[i][j][k] *= 2.0 / (1.0 - float64(j)/float64(sf.P) - float64(k)/float64(sf.P) + 1e-10)
				case 1: // s-derivative
					DU[i][j][k] *= 2.0 / (1.0 - float64(k)/float64(sf.P) + 1e-10)
				case 2: // t-derivative
					DU[i][j][k] *= 2.0
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

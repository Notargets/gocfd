package DG2D

import (
	"fmt"
	"math"

	"github.com/notargets/gocfd/utils"
	"gonum.org/v1/gonum/mat"
)

// monomialVandermonde builds an n×M Vandermonde matrix (using a monomial basis)
// from the provided nodes. For each node x (assumed to be in [-1,1]),
// it computes: 1, x, x^2, …, x^(M-1).
func monomialVandermonde(nodes []float64, M int) utils.Matrix {
	n := len(nodes)
	V := utils.NewMatrix(n, M)
	for i, x := range nodes {
		for j := 0; j < M; j++ {
			V.Set(i, j, math.Pow(x, float64(j)))
		}
	}
	return V
}

// EdgeInterpolator precomputes the Vandermonde (and its inverse)
// for a fixed set of sample nodes (all assumed in [-1,1]).
type EdgeInterpolator struct {
	M                int          // number of sample nodes (and thus polynomial degree = M-1)
	VSample          utils.Matrix // Vandermonde matrix (M×M) built from sample nodes.
	VSampleInv       utils.Matrix // Inverse of the Vandermonde matrix.
	lambdaC, lambdaH float64      // Exponential Decay scaling for const and HOT
	Nodes            []float64    // The sample node locations.
}

// NewEdgeInterpolator constructs an EdgeInterpolator using the given sampleNodes.
// Here M is taken to be len(sampleNodes).
// func NewEdgeInterpolator(sampleNodes []float64, LambdaC, LambdaH float64) *EdgeInterpolator {
func NewEdgeInterpolator(sampleNodes []float64, consts ...float64) *EdgeInterpolator {
	var (
		lambdaC, lambdaH float64
	)
	lambdaC = 2.0 // Controls exponential decay of constant term.
	lambdaH = 2.0 // Controls exponential decay of higher order terms.
	switch len(consts) {
	case 0:
		// keep additional default constants
	case 1:
		lambdaC = consts[0] // Controls exponential decay of constant term.
	case 2:
		lambdaC = consts[0] // Controls exponential decay of constant term.
		lambdaH = consts[1] // Controls exponential decay of higher order terms.
	}
	M := len(sampleNodes)
	V := monomialVandermonde(sampleNodes, M)
	VInv := V.InverseWithCheck() // Inverse() returns a pointer to a new matrix.
	return &EdgeInterpolator{
		M:          M,
		VSample:    V,
		VSampleInv: VInv,
		Nodes:      sampleNodes,
		lambdaC:    lambdaC,
		lambdaH:    lambdaH,
	}
}

// Polynomial represents a polynomial in the monomial basis:
// P(x) = a₀ + a₁ x + a₂ x² + … + aₘ₋₁ x^(M-1).
type Polynomial struct {
	Coeffs []float64 // Coefficients [a₀, a₁, …, aₘ₋₁]
	M      int       // Number of coefficients.
}

// Evaluate returns P(x) = Σ₀^(M-1) a_j * x^j.
func (p *Polynomial) Evaluate(x float64) float64 {
	sum := 0.0
	xpow := 1.0
	for j := 0; j < p.M; j++ {
		sum += p.Coeffs[j] * xpow
		xpow *= x
	}
	return sum
}

// Derivative returns P'(x) = Σ₁^(M-1) j * a_j * x^(j-1).
func (p *Polynomial) Derivative(x float64) float64 {
	sum := 0.0
	xpow := 1.0
	for j := 1; j < p.M; j++ {
		sum += float64(j) * p.Coeffs[j] * xpow
		xpow *= x
	}
	return sum
}

// modulatePolynomialEx applies two exponential filters:
//   - The constant term is moved toward mid as:
//     a₀′ = mid + (a₀ – mid) * exp(-lambdaC * (1 - fScale))
//   - The higher–order terms (j>=1) are filtered as:
//     aⱼ′ = aⱼ * exp(-lambdaH * (1 - fScale))
func modulatePolynomialEx(poly *Polynomial, fScale, mid, lambdaC, lambdaH float64) *Polynomial {
	aMod := make([]float64, poly.M)
	aMod[0] = mid + (poly.Coeffs[0]-mid)*math.Exp(-lambdaC*(1-fScale))
	for j := 1; j < poly.M; j++ {
		aMod[j] = poly.Coeffs[j] * math.Exp(-lambdaH*(1-fScale))
	}
	return &Polynomial{
		Coeffs: aMod,
		M:      poly.M,
	}
}

// analyticalExtrema computes the minimum and maximum values of poly(x) over x in [-1,1].
// It does so by evaluating the polynomial at the endpoints and at the critical points (roots of P'(x)=0).
func analyticalExtrema(poly *Polynomial) (float64, float64) {
	candidates := []float64{-1, 1}
	roots := findCriticalPointsRobust(poly)
	candidates = append(candidates, roots...)
	minVal := poly.Evaluate(candidates[0])
	maxVal := minVal
	for _, x := range candidates {
		val := poly.Evaluate(x)
		if val < minVal {
			minVal = val
		}
		if val > maxVal {
			maxVal = val
		}
	}
	return minVal, maxVal
}

// FitAndBoundPolynomial fits a polynomial (via the precomputed Vandermonde inverse)
// from the provided sampleValues (length M) and then modulates the coefficients so that
// the polynomial'S value over [-1,1] remains within [VMin, VMax].
// If the initial polynomial already meets the bounds, a message is printed and the original
// polynomial is returned.
// Otherwise, a binary search on fScale (in [0,1]) is performed, using an exponential modulation
// on both the constant and higher-order coefficients.
func (ei *EdgeInterpolator) FitAndBoundPolynomial(sampleValues []float64, VMin, VMax float64) *Polynomial {
	if len(sampleValues) != ei.M {
		panic("length of sampleValues must equal number of sample nodes")
	}
	// Create the right-hand side vector f.
	f := utils.NewMatrix(ei.M, 1)
	for i, val := range sampleValues {
		f.Set(i, 0, val)
	}
	// Compute coefficients: a = VSampleInv * f.
	aMat := ei.VSampleInv.Mul(f) // M×1 result.
	a := make([]float64, ei.M)
	for i := 0; i < ei.M; i++ {
		a[i] = aMat.At(i, 0)
	}
	// Build the initial polynomial.
	poly := &Polynomial{
		Coeffs: a,
		M:      ei.M,
	}
	// Compute the midpoint.
	mid := (VMin + VMax) / 2

	// Check initial extrema.
	pMin, pMax := analyticalExtrema(poly)
	if pMin >= VMin && pMax <= VMax {
		fmt.Printf("Initial polynomial is within bounds: min=%f, max=%f\n", pMin, pMax)
		return poly
	}

	// Binary search for fScale in [0,1].
	tol := 1e-6
	low, high := 0.0, 1.0
	var midScale float64
	var modPoly *Polynomial
	for high-low > tol {
		midScale = (low + high) / 2.0
		modPoly = modulatePolynomialEx(poly, midScale, mid, ei.lambdaC, ei.lambdaH)
		modMin, modMax := analyticalExtrema(modPoly)
		if modMin >= VMin && modMax <= VMax {
			// The modulated polynomial meets the bounds; try allowing more high-order content.
			low = midScale
		} else {
			high = midScale
		}
	}
	return modulatePolynomialEx(poly, low, mid, ei.lambdaC, ei.lambdaH)
}

// findCriticalPointsRobust returns the real roots of P'(x)=0 in [-1,1].
// For M <= 3 (polynomials up to quadratic) closed–form formulas are used.
// For 4 <= M <= 7 a companion–matrix approach is applied.
func findCriticalPointsRobust(poly *Polynomial) []float64 {
	if poly.M <= 2 {
		// Linear polynomial: derivative is constant.
		return []float64{}
	} else if poly.M == 3 {
		// Quadratic polynomial: P(x)= a0 + a1 x + a2 x^2, derivative: a1 + 2*a2 x.
		a1 := poly.Coeffs[1]
		a2 := poly.Coeffs[2]
		if math.Abs(a2) < 1e-12 {
			return []float64{}
		}
		root := -a1 / (2 * a2)
		if root >= -1 && root <= 1 {
			return []float64{root}
		}
		return []float64{}
	} else {
		// For poly.M between 4 and 7, the derivative P'(x) is a polynomial of degree M-2.
		// Build its coefficients in the monomial basis:
		// Let D(x) = Σ_{k=0}^{d} c[k] x^k, with d = poly.M-2 and c[k] = (k+1)*a[k+1].
		d := poly.M - 2
		c := make([]float64, d+1)
		for k := 0; k <= d; k++ {
			c[k] = float64(k+1) * poly.Coeffs[k+1]
		}
		// If the leading coefficient is nearly zero, return no roots.
		if math.Abs(c[d]) < 1e-12 {
			return []float64{}
		}
		// Normalize to make the polynomial monic.
		coeffs := make([]float64, d)
		for i := 0; i < d; i++ {
			coeffs[i] = c[i] / c[d]
		}
		// Construct the companion matrix of size d×d.
		comp := mat.NewDense(d, d, nil)
		// The first row is [-coeffs[d-1], -coeffs[d-2], …, -coeffs[0]].
		for j := 0; j < d; j++ {
			comp.Set(0, j, -coeffs[d-1-j])
		}
		// SetScalar the subdiagonal to 1.
		for i := 1; i < d; i++ {
			comp.Set(i, i-1, 1.0)
		}
		// Compute the eigenvalues of the companion matrix.
		var eig mat.Eigen
		if ok := eig.Factorize(comp, mat.EigenRight); !ok {
			return []float64{}
		}
		vals := eig.Values(nil)
		var roots []float64
		for _, lambda := range vals {
			if math.Abs(imag(lambda)) < 1e-8 {
				r := real(lambda)
				if r >= -1 && r <= 1 {
					roots = append(roots, r)
				}
			}
		}
		return roots
	}
}

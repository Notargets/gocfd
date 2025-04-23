package DG2D

import (
	"fmt"
	"math"
	"math/rand"
	"sort"
	"time"

	"gonum.org/v1/gonum/optimize"

	"github.com/notargets/gocfd/utils"

	"gonum.org/v1/gonum/mat"
)

// RS holds a point in the standard (r,s) triangle
// with vertices (-1,-1),(1,-1),(-1,1).
type Point struct{ R, S float64 }

// monomialExponents returns all (i,j) with i+j ≤ deg.
func monomialExponents(deg int) [][2]int {
	exps := make([][2]int, 0, (deg+1)*(deg+2)/2)
	for total := 0; total <= deg; total++ {
		for i := 0; i <= total; i++ {
			exps = append(exps, [2]int{i, total - i})
		}
	}
	return exps
}

// uniformBarycentric returns barycentric triples on (0,0)-(1,0)-(0,1)
// at resolution n: (i/n, j/n, 1−i/n−j/n).
func uniformBarycentric(n int) [][3]float64 {
	out := make([][3]float64, 0, (n+1)*(n+2)/2)
	inv := 1.0 / float64(n)
	for i := 0; i <= n; i++ {
		li := float64(i) * inv
		for j := 0; j <= n-i; j++ {
			lj := float64(j) * inv
			lk := 1 - li - lj
			out = append(out, [3]float64{li, lj, lk})
		}
	}
	return out
}

// logDet returns ln|det(A)|.  It panics only if the determinant is zero.
func logDet(A utils.Matrix) float64 {
	var lu mat.LU
	lu.Factorize(A)
	logAbs, sign := lu.LogDet()
	if sign == 0 {
		// exact zero determinant
		panic("logDet: zero determinant")
	}
	// logAbs may be negative—that's okay
	return logAbs
}

// UniformRS returns a “barycentric” grid of resolution n
// mapped straight into the (r,s) triangle.  You get exactly
// (n+1)(n+2)/2 points whose little cells all have equal area.
func UniformRS(n int) []Point {
	if n < 1 {
		return nil
	}
	pts := make([]Point, 0, (n+1)*(n+2)/2)
	inv := 1.0 / float64(n)
	for i := 0; i <= n; i++ {
		li := float64(i) * inv // barycentric λ₁
		for j := 0; j <= n-i; j++ {
			lj := float64(j) * inv // λ₂
			lk := 1 - li - lj      // λ₃
			// map λ₂,λ₃ in [0,1] to r,s in [-1,1]:
			//    r = 2*λ₂ - 1,  s = 2*λ₃ - 1
			pts = append(pts, Point{
				R: 2*lj - 1,
				S: 2*lk - 1,
			})
		}
	}
	return pts
}

func getInteriorPoints(pts []Point) (intPts []Point) {
	// var interiorCount, extCount int
	// var vertexCount int
	near := func(x, y float64) bool {
		var (
			tol = 1.e-6
		)
		if math.Abs(x-y) < tol {
			return true
		}
		return false
	}
	var isExterior bool
	for _, pt := range pts {
		isExterior = false
		rr, ss := pt.R, pt.S
		if near(ss, -1) || near(rr, -1) ||
			near(ss, 1) || near(rr, 1) ||
			near(rr+ss, 0) {
			// extCount++
			isExterior = true
		}
		if (near(ss, -1) && near(rr, -1)) ||
			(near(ss, 1) && near(rr, -1)) ||
			(near(ss, -1) && near(rr, 1)) {
			// vertexCount++
			isExterior = true
		}
		if !isExterior {
			intPts = append(intPts, pt)
		}
	}
	// totalCount := len(pts)
	// interiorCount = totalCount - extCount
	// singleEdge := (totalCount - interiorCount - vertexCount) / 3
	// _ = singleEdge
	return
}

// OptimizeRTNodesCont does a continuous off‐lattice Fekete‐style interior
// optimization for the RT element of order P, with Gauss–Legendre edges.
func OptimizeRTNodesCont(P int) ([]Point, error) {
	if P < 1 {
		return nil, fmt.Errorf("P must be ≥1, got %d", P)
	}

	// 1) Compute the P+1 Gauss–Legendre nodes on [-1,1].
	gl := gaussLegendreNodes(P + 1)
	// Drop the two endpoints (-1 and +1) to get exactly P interior-edge points.
	// tvals := gl[1 : len(gl)-1] // len(tvals) == P

	// 2) Lay them out on the three edges: bottom, left, hypotenuse.
	// Drop none — Legendre roots already lie strictly in (-1,1).
	// Place them on bottom (s=-1), left (r=-1), and hypotenuse (s=-r).
	edges := make([]Point, 0, 3*len(gl))
	for _, t := range gl {
		// bottom edge:  (r ∈ (-1,1), s = -1)
		edges = append(edges, Point{R: t, S: -1})
		// left edge:    (r = -1, s ∈ (-1,1))
		edges = append(edges, Point{R: -1, S: t})
		// hypotenuse:   line s = -r
		edges = append(edges, Point{R: t, S: -t})
	}

	// 3) Prepare monomial exponents for total degree ≤ P-1.
	M := P * (P + 1) / 2
	exps := monomialExponents(P - 1)
	if len(exps) != M {
		return nil, fmt.Errorf("internal: expected %d monomials, got %d", M, len(exps))
	}

	// 4) Define the optimization problem over 2*M real parameters.
	prob := optimize.Problem{
		Func: func(x []float64) float64 {
			// Build the M×M Vandermonde for the interior points.
			PhiI := utils.NewMatrix(M, M)
			for i := 0; i < M; i++ {
				xi, yi := x[2*i], x[2*i+1]
				l1 := logistic(xi)
				l2 := logistic(yi) * (1 - l1)
				l3 := 1 - l1 - l2
				r, s := 2*l2-1, 2*l3-1
				for j, e := range exps {
					PhiI.Set(i, j, math.Pow(r, float64(e[0]))*math.Pow(s, float64(e[1])))
				}
			}
			// Gram = ΦᵀΦ
			Gram := PhiI.Transpose().Mul(PhiI)
			// LU log‐det
			// 3) LU‐logdet
			var lu mat.LU
			lu.Factorize(Gram)
			logAbs, sign := lu.LogDet()
			if sign == 0 {
				return math.Inf(1)
			}

			// 4) hinge‐penalty for boundary proximity (as before)
			const (
				eps   = 1e-1 // threshold
				alpha = 1e3  // strength
				// gamma controls how strongly you push for conditioning
				gamma = 1e-0
			)
			penalty := 0.0
			for i := 0; i < M; i++ {
				xi, yi := x[2*i], x[2*i+1]
				l1 := logistic(xi)
				l2 := logistic(yi) * (1 - l1)
				l3 := 1 - l1 - l2

				if l1 < eps {
					d := eps - l1
					penalty += alpha * d * d
				}
				if l2 < eps {
					d := eps - l2
					penalty += alpha * d * d
				}
				if l3 < eps {
					d := eps - l3
					penalty += alpha * d * d
				}
			}

			// 5) condition‐number penalty: use general Eigen on the dense Gram
			var eig mat.Eigen
			var GramDense mat.Dense
			GramDense.CloneFrom(Gram) // copy your utils.Matrix into a *mat.Dense

			// Factorize asking for NO eigenvectors
			ok := eig.Factorize(&GramDense, mat.EigenNone)
			if !ok {
				return math.Inf(1) // something went very wrong
			}

			// Grab the (complex128) eigenvalues, extract reals, sort, and build κ
			raw := eig.Values(nil) // []complex128
			reals := make([]float64, len(raw))
			for i, c := range raw {
				reals[i] = real(c)
			}
			sort.Float64s(reals) // now sorted ascending
			lambdaMin, lambdaMax := reals[0], reals[len(reals)-1]

			logCond := math.Log(lambdaMax / lambdaMin)
			return -logAbs + penalty + gamma*logCond

		},
	}

	// 5) Initialize x with small random jitter so we start nonsingular.
	initX := make([]float64, 2*M)
	rand.Seed(time.Now().UnixNano())
	for i := range initX {
		initX[i] = rand.Float64()*2 - 1
	}

	// 6) Run Nelder–Mead with default convergence.
	settings := &optimize.Settings{Recorder: nil}
	result, err := optimize.Minimize(prob, initX, settings, &optimize.NelderMead{})
	if err != nil {
		return nil, fmt.Errorf("optimization failed: %v", err)
	}

	// 7) Reconstruct the optimized interior (r,s) points.
	xOpt := result.X
	interior := make([]Point, 0, M)
	for i := 0; i < M; i++ {
		l1 := logistic(xOpt[2*i])
		l2 := logistic(xOpt[2*i+1]) * (1 - l1)
		l3 := 1 - l1 - l2
		interior = append(interior, Point{R: 2*l2 - 1, S: 2*l3 - 1})
	}

	// 8) Return the fixed edges first, then the optimized interior.
	return append(edges, interior...), nil
}

// logistic maps ℝ→(0,1).
func logistic(x float64) float64 {
	return 1 / (1 + math.Exp(-x))
}

// gaussLegendreNodes returns the n roots of Pₙ(x) on [-1,1].
func gaussLegendreNodes(n int) []float64 {
	if n == 1 {
		return []float64{0}
	}
	J := mat.NewSymDense(n, nil)
	for i := 0; i < n-1; i++ {
		k := float64(i + 1)
		v := k / math.Sqrt(4*k*k-1)
		J.SetSym(i, i+1, v)
		J.SetSym(i+1, i, v)
	}
	var eig mat.EigenSym
	eig.Factorize(J, true)
	nodes := eig.Values(nil)
	sort.Float64s(nodes)
	return nodes
}

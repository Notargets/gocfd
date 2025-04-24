package DG2D

import (
	"fmt"
	"math"
	"math/rand"
	"sort"
	"strconv"

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
func MakeEdgePointsFromDist(edgeDist []float64) (pts []Point) {
	var (
		le = len(edgeDist)
	)
	pts = make([]Point, 3*le)
	for i, r := range edgeDist {
		pts[i] = Point{R: r, S: -1}
		ii := i + le
		pts[ii] = Point{R: -r, S: r}
		iii := i + 2*le
		pts[iii] = Point{R: -1, S: -r}
	}
	return
}

func AddEdgePointsToInterior(edgeDist []float64, R, S utils.Vector) (pts []Point) {
	interiorPoints := make([]Point, R.Len())
	for i, r := range R.DataP {
		interiorPoints[i] = Point{R: r, S: S.DataP[i]}
	}
	edgePoints := MakeEdgePointsFromDist(edgeDist)
	pts = append(edgePoints, interiorPoints...)
	return
}

func MakeRSFromPoints(pts []Point) (R, S utils.Vector) {
	R = utils.NewVector(len(pts))
	S = utils.NewVector(len(pts))
	for i, pt := range pts {
		R.DataP[i] = pt.R
		S.DataP[i] = pt.S
	}
	return
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

// UniformRSAlpha returns a “barycentric” grid of resolution n
// mapped into the (r,s) triangle and then uniformly shrunk
// toward the triangle’s centroid by factor α∈[0,1].
// You get exactly (n+1)(n+2)/2 points; α=1 ⇒ no shrink; α=0 ⇒ all at centroid.
func UniformRSAlpha(n int, alpha float64) []Point {
	// reference‐triangle centroid in (r,s):
	const rc, sc = -1.0 / 3.0, -1.0 / 3.0
	if n == 0 {
		return []Point{Point{rc, sc}}
	} else if n < 0 {
		panic("n must be positive, got " + strconv.Itoa(n))
	}

	pts := make([]Point, 0, (n+1)*(n+2)/2)
	inv := 1.0 / float64(n)

	for i := 0; i <= n; i++ {
		li := float64(i) * inv // λ₁
		for j := 0; j <= n-i; j++ {
			lj := float64(j) * inv // λ₂
			lk := 1 - li - lj      // λ₃

			// map (λ₂,λ₃) → (r,s)
			r := 2*lj - 1
			s := 2*lk - 1

			// shrink toward centroid (rc,sc)
			r = rc + alpha*(r-rc)
			s = sc + alpha*(s-sc)

			pts = append(pts, Point{R: r, S: s})
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

// OptimizeRTWithParams does off‐lattice interior optimization
// for RT‐P in the standard triangle, with fixed GL edges.
// eps   = minimum barycentric threshold,
// alpha = hinge strength,
// gamma = conditioning weight.
func OptimizeRTWithParams(
	P int,
	eps, alpha, gamma float64,
) ([]Point, error) {
	// Objective J(x):
	//
	//   J(x) = – log │det(Φ(x)ᵀ · Φ(x))│
	//          + α ∑_{i=1}^M ∑_{j=1}^3 (ε – λᵢⱼ)_+²
	//          + γ · log(λ_max / λ_min)
	//
	// where
	//
	//   • Φ(x)    is the M×M Vandermonde matrix built from the M interior points.
	//   • λᵢ₁, λᵢ₂, λᵢ₃  are the three barycentric coordinates (ℓ₁,ℓ₂,ℓ₃) of point i.
	//   • (y)_+  = max(0, y)  is the hinge function.
	//   • λ_min, λ_max  are the smallest and largest eigenvalues of the Gram matrix ΦᵀΦ.
	//   • ε  is the minimum‐barycentric threshold.
	//   • α  is the hinge strength on barycentric bounds.
	//   • γ  is the weight on the log‐condition‐number penalty.
	if P < 1 {
		return nil, fmt.Errorf("P must be ≥1, got %d", P)
	}

	// --- build fixed GL edge nodes (these never move) ---
	gl := gaussLegendreNodes(P + 1)
	edges := make([]Point, 0, 3*len(gl))
	for _, t := range gl {
		edges = append(edges,
			Point{R: t, S: -1}, // bottom
			Point{R: 1, S: t},  // right
			Point{R: -1, S: t}, // left
		)
	}

	// interior count and monomial exponents
	M := P * (P + 1) / 2
	exps := monomialExponents(P - 1)

	// problem definition
	prob := optimize.Problem{
		Func: func(x []float64) float64 {
			// 1) build Vandermonde Φ (M×M) for current interior guess x
			Phi := utils.NewMatrix(M, M)
			for i := 0; i < M; i++ {
				xi, yi := x[2*i], x[2*i+1]
				// stick-breaking → barycentric
				l1 := logistic(xi)
				l2 := logistic(yi) * (1.0 - l1)
				l3 := 1.0 - l1 - l2
				// map to (r,s)
				r, s := 2*l2-1, 2*l3-1

				for j, e := range exps {
					Phi.Set(i, j,
						math.Pow(r, float64(e[0]))*
							math.Pow(s, float64(e[1])),
					)
				}
			}

			// 2) Gram = ΦᵀΦ
			Gram := Phi.Transpose().Mul(Phi)

			// 3) LU‐logdet part
			var lu mat.LU
			lu.Factorize(Gram)
			logAbs, sign := lu.LogDet()
			if sign == 0 {
				return math.Inf(1) // singular → infinite penalty
			}
			// core objective
			obj := -logAbs

			// 4) boundary‐hinge penalty on each barycentric λ₁,₂,₃
			for i := 0; i < M; i++ {
				xi, yi := x[2*i], x[2*i+1]
				l1 := logistic(xi)
				l2 := logistic(yi) * (1.0 - l1)
				l3 := 1.0 - l1 - l2

				for _, l := range []float64{l1, l2, l3} {
					if l < eps {
						d := eps - l
						obj += alpha * d * d
					}
				}
			}

			// 5) conditioning penalty via Gram eigenvalues
			// copy into a *mat.Dense to call Eigen
			var Gmat mat.Dense
			Gmat.CloneFrom(Gram)
			var eig mat.Eigen
			ok := eig.Factorize(&Gmat, mat.EigenNone)
			if !ok {
				return math.Inf(1)
			}
			vals := eig.Values(nil)
			reals := make([]float64, len(vals))
			for i, c := range vals {
				reals[i] = real(c)
			}
			sort.Float64s(reals)
			λmin, λmax := reals[0], reals[len(reals)-1]
			// log‐condition
			obj += gamma * math.Log(λmax/λmin)

			return obj
		},
	}

	// 6) initialize x with small random in [-1,1]
	initX := make([]float64, 2*M)
	for i := range initX {
		initX[i] = rand.Float64()*2 - 1
	}

	// 7) run Nelder–Mead
	settings := &optimize.Settings{Recorder: nil}
	result, err := optimize.Minimize(prob, initX, settings, &optimize.NelderMead{})
	if err != nil {
		return nil, fmt.Errorf("optimize failed: %v", err)
	}

	// 8) unpack final interior points
	Xopt := result.X
	interior := make([]Point, 0, M)
	for i := 0; i < M; i++ {
		l1 := logistic(Xopt[2*i])
		l2 := logistic(Xopt[2*i+1]) * (1.0 - l1)
		l3 := 1.0 - l1 - l2
		interior = append(interior, Point{
			R: 2*l2 - 1,
			S: 2*l3 - 1,
		})
	}

	// 9) return edges followed by interior
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

// AnalyzeRTOutput takes the full output of OptimizeRTNodesCont:
//
//	pts = [ e0,e1,…,e_{3(P+1)-1},  i0,i1,…,i_{M-1} ]
//
// and returns:
//
//	cond    = condition number of the M×M Vandermonde at the interior points,
//	lebEdge = Lebesgue constant *on the fixed edge nodes* e_j.
func AnalyzeRTOutput(P int, pts []Point) (cond, lebEdge float64) {
	E := 3 * (P + 1)
	M := P * (P + 1) / 2

	if len(pts) != E+M {
		panic(fmt.Errorf(
			"AnalyzeRTOutput: expected %d total points (3*(P+1)+P(P+1)/2), got %d",
			E+M, len(pts),
		))
	}

	// split edges vs interior
	edges := pts[:E]
	interior := pts[E:]

	// 1) build monomial exponents for total deg ≤ P-1
	exps := monomialExponents(P - 1)

	// 2) build the M×M Vandermonde V for the interior points
	V := mat.NewDense(M, M, nil)
	for i, p := range interior {
		for j, e := range exps {
			V.Set(i, j,
				math.Pow(p.R, float64(e[0]))*
					math.Pow(p.S, float64(e[1])),
			)
		}
	}

	// 3) condition number via thin SVD
	var svd mat.SVD
	ok := svd.Factorize(V, mat.SVDThin)
	if !ok {
		panic(fmt.Errorf("AnalyzeRTOutput: SVD failed on Vandermonde"))
	}
	sig := svd.Values(nil)
	sort.Float64s(sig)
	cond = sig[len(sig)-1] / sig[0]

	// 4) LU‐factor Vandermonde for fast solves
	var lu mat.LU
	lu.Factorize(V)

	// 5) compute edge‐Lebesgue constant:
	//    for each fixed edge node e in edges, build φ(e), solve V α = φ,
	//    sum |α_i|, and take the max.
	lebEdge = 0
	for _, ept := range edges {
		// build φ(ept)
		phi := mat.NewVecDense(M, nil)
		for j, ex := range exps {
			phi.SetVec(j,
				math.Pow(ept.R, float64(ex[0]))*
					math.Pow(ept.S, float64(ex[1])),
			)
		}
		// solve for α: V·α = φ
		var alpha mat.VecDense
		alpha.SolveVec(&lu, phi)

		sum := 0.0
		for i := 0; i < M; i++ {
			sum += math.Abs(alpha.AtVec(i))
		}
		if sum > lebEdge {
			lebEdge = sum
		}
	}

	return cond, lebEdge
}

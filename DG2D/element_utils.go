package DG2D

import (
	"math"

	"github.com/notargets/gocfd/DG1D"
	"github.com/notargets/gocfd/utils"
)

// Purpose  : Compute (x,y) nodes in equilateral triangle for
//
//	polynomial of order N
func Nodes2D(N int) (x utils.Vector, y utils.Vector) {
	var (
		alpha                                                               float64
		Np                                                                  = (N + 1) * (N + 2) / 2
		L1, L2, L3                                                          utils.Vector
		blend1, blend2, blend3, warp1, warp2, warp3, warpf1, warpf2, warpf3 []float64
	)
	L1, L2, L3, x, y =
		utils.NewVector(Np), utils.NewVector(Np), utils.NewVector(Np), utils.NewVector(Np), utils.NewVector(Np)
	l1d, l2d, l3d, xd, yd := L1.DataP, L2.DataP, L3.DataP, x.DataP, y.DataP
	blend1, blend2, blend3, warp1, warp2, warp3 =
		make([]float64, Np), make([]float64, Np), make([]float64, Np), make([]float64, Np), make([]float64, Np), make([]float64, Np)

	alpopt := []float64{
		0.0000, 0.0000, 1.4152, 0.1001, 0.2751,
		0.9800, 1.0999, 1.2832, 1.3648, 1.4773,
		1.4959, 1.5743, 1.5770, 1.6223, 1.6258,
	}
	if N < 16 {
		alpha = alpopt[N-1]
	} else {
		alpha = 5. / 3.
	}
	// Create equidistributed nodes on equilateral triangle
	fn := 1. / float64(N)
	var sk int
	for n := 0; n < N+1; n++ {
		for m := 0; m < (N + 1 - n); m++ {
			l1d[sk] = float64(n) * fn
			l3d[sk] = float64(m) * fn
			sk++
		}
	}
	for i := range xd {
		l2d[i] = 1 - l1d[i] - l3d[i]
		xd[i] = l3d[i] - l2d[i]
		yd[i] = (2*l1d[i] - l3d[i] - l2d[i]) / math.Sqrt(3)
		// Compute blending function at each node for each edge
		blend1[i] = 4 * l2d[i] * l3d[i]
		blend2[i] = 4 * l1d[i] * l3d[i]
		blend3[i] = 4 * l1d[i] * l2d[i]
	}
	// Amount of warp for each node, for each edge
	warpf1 = Warpfactor(N, L3.Copy().Subtract(L2))
	warpf2 = Warpfactor(N, L1.Copy().Subtract(L3))
	warpf3 = Warpfactor(N, L2.Copy().Subtract(L1))
	// Combine blend & warp
	for i := range warpf1 {
		warp1[i] = blend1[i] * warpf1[i] * (1 + utils.POW(alpha*l1d[i], 2))
		warp2[i] = blend2[i] * warpf2[i] * (1 + utils.POW(alpha*l2d[i], 2))
		warp3[i] = blend3[i] * warpf3[i] * (1 + utils.POW(alpha*l3d[i], 2))
	}
	// Accumulate deformations associated with each edge
	for i := range xd {
		xd[i] += warp1[i] + math.Cos(2*math.Pi/3)*warp2[i] + math.Cos(4*math.Pi/3)*warp3[i]
		yd[i] += math.Sin(2*math.Pi/3)*warp2[i] + math.Sin(4*math.Pi/3)*warp3[i]
	}
	return
}

func Warpfactor(N int, rout utils.Vector) (warpF []float64) {
	var (
		Nr   = rout.Len()
		Pmat = utils.NewMatrix(N+1, Nr)
	)
	// Compute LGL and equidistant node distribution
	LGLr := DG1D.JacobiGL(0, 0, N)
	req := utils.NewVector(N+1).Linspace(-1, 1)
	Veq := DG1D.Vandermonde1D(N, req)
	// Evaluate Lagrange polynomial at rout
	for i := 0; i < (N + 1); i++ {
		Pmat.M.SetRow(i, DG1D.JacobiP(rout, 0, 0, i))
	}
	Lmat := Veq.Transpose().LUSolve(Pmat)
	// Compute warp factor
	warp := Lmat.Transpose().Mul(LGLr.Subtract(req).ToMatrix())
	// Scale factor
	zerof := rout.Copy().Apply(func(val float64) (res float64) {
		if math.Abs(val) < (1.0 - (1e-10)) {
			res = 1.
		}
		return
	})
	sf := zerof.Copy().ElMul(rout).Apply(func(val float64) (res float64) {
		res = 1 - val*val
		return
	})
	w2 := warp.Copy()
	warp.ElDiv(sf.ToMatrix()).Add(w2.ElMul(zerof.AddScalar(-1).ToMatrix()))
	warpF = warp.DataP
	return
}

func RStoAB(R, S utils.Vector) (a, b utils.Vector) {
	var (
		Np     = R.Len()
		rd, sd = R.DataP, S.DataP
	)
	ad, bd := make([]float64, Np), make([]float64, Np)
	for n, sval := range sd {
		/*
			if sval != 1 {
				ad[n] = 2*(1+rd[n])/(1-sval) - 1
			} else {
				ad[n] = -1
			}
			bd[n] = sval
		*/
		ad[n], bd[n] = rsToab(rd[n], sval)
	}
	a, b = utils.NewVector(Np, ad), utils.NewVector(Np, bd)
	return
}

func rsToab(r, s float64) (a, b float64) {
	if s != 1 {
		a = 2*(1+r)/(1-s) - 1
	} else {
		a = -1
	}
	b = s
	return
}

// function [R,S] = xytors(x,y)
// Purpose : Transfer from (x,y) in equilateral triangle
//
//	to (R,S) coordinates in standard triangle
func XYtoRS(x, y utils.Vector) (r, s utils.Vector) {
	r, s = utils.NewVector(x.Len()), utils.NewVector(x.Len())
	var (
		xd, yd = x.DataP, y.DataP
		rd, sd = r.DataP, s.DataP
	)
	sr3 := math.Sqrt(3)
	for i := range xd {
		l1 := (sr3*yd[i] + 1) / 3
		l2 := (-3*xd[i] - sr3*yd[i] + 2) / 6
		l3 := (3*xd[i] - sr3*yd[i] + 2) / 6
		rd[i] = -l2 + l3 - l1
		sd[i] = -l2 - l3 + l1
	}
	return
}

func CalculateElementLocalGeometry(EToV utils.Matrix, VX, VY, R, S utils.Vector) (X, Y utils.Matrix) {
	/*
		For input values of vector field [R,S], transform them into element local [X,Y]
	*/
	va, vb, vc := EToV.Col(0), EToV.Col(1), EToV.Col(2)
	X = R.Copy().Add(S).Scale(-1).Outer(VX.SubsetIndex(va.ToIndex())).Add(
		R.Copy().AddScalar(1).Outer(VX.SubsetIndex(vb.ToIndex()))).Add(
		S.Copy().AddScalar(1).Outer(VX.SubsetIndex(vc.ToIndex()))).Scale(0.5)
	Y = R.Copy().Add(S).Scale(-1).Outer(VY.SubsetIndex(va.ToIndex())).Add(
		R.Copy().AddScalar(1).Outer(VY.SubsetIndex(vb.ToIndex()))).Add(
		S.Copy().AddScalar(1).Outer(VY.SubsetIndex(vc.ToIndex()))).Scale(0.5)
	return
}

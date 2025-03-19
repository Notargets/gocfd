package DG2D

import (
	"math"

	"github.com/notargets/gocfd/utils"
)

type ModeAliasShockFinder struct {
	Element *LagrangeElement2D
	Clipper utils.Matrix // Matrix used to clip the topmost mode from the solution polynomial, used in shockfinder
	Np      int
	Q, Qalt utils.Matrix // scratch storage for evaluating the moment
	Kappa   float64
}

func (dfr *DFR2D) NewAliasShockFinder(Kappa float64) (sf *ModeAliasShockFinder) {
	var (
		element = dfr.SolutionElement
		Np      = element.Np
		N       = element.N
	)
	sf = &ModeAliasShockFinder{
		Element: element,
		Np:      Np,
		Q:       utils.NewMatrix(Np, 1),
		Qalt:    utils.NewMatrix(Np, 1),
		Kappa:   Kappa,
	}
	/*
		The "Clipper" matrix drops the last mode from the polynomial and forms an alternative field of values at the node
		points based on a polynomial with one less term. In other words, if we have a polynomial of degree "p", expressed
		as values at Np node points, multiplying the Node point values vector by Clipper produces an alternative version
		of the node values based on truncating the last polynomial mode.
	*/
	sf.Clipper =
		element.JB2D.V.Mul(dfr.CutoffFilter2D(N, N, 0)).Mul(element.JB2D.Vinv)
	return
}

func (dfr *DFR2D) CutoffFilter2D(N, NCutoff int, frac float64) (diag utils.Matrix) {
	/*
		The NCutoff is inclusive, so if you want to clip the top mode at N, input N
	*/
	var (
		Np = (N + 1) * (N + 2) / 2
	)
	data := make([]float64, Np)
	for ii := 0; ii < Np; ii++ {
		data[ii] = 1.
	}
	var ii int
	for i := 0; i <= N; i++ {
		for j := 0; j <= N-i; j++ {
			if i+j >= NCutoff {
				data[ii] = frac
			}
			ii++
		}
	}
	diag = utils.NewDiagMatrix(Np, data)
	return
}

func (sf *ModeAliasShockFinder) ElementHasShock(q []float64) (i bool) {
	if sf.ShockIndicator(q) > 0.0075 {
		i = true
	}
	return
}

func (sf *ModeAliasShockFinder) ShockIndicator(q []float64) (sigma float64) {
	/*
		Original method by Persson, constants chosen to match Zhiqiang, et. al.
	*/
	var (
		Se          = math.Log10(sf.moment(q))
		kappa       = sf.Kappa
		S0          = 4. / math.Pow(float64(sf.Element.N), 4.)
		left, right = S0 - kappa, S0 + kappa
		ookappa     = 0.5 / kappa
	)
	switch {
	case Se < left:
		sigma = 0.
	case Se >= left && Se <= right:
		sigma = 0.5 * (1. + math.Sin(math.Pi*ookappa*(Se-S0)))
	case Se > right:
		sigma = 1.
	}
	return
}

func (sf *ModeAliasShockFinder) moment(q []float64) (m float64) {
	var (
		Np            = sf.Np
		U, UClipped   = sf.Q, sf.Qalt
		UD, UClippedD = U.DataP, UClipped.DataP
		MD            = sf.Element.MassMatrix.DataP
	)
	copy(sf.Q.DataP, q)
	/*
		Evaluate the L2 moment of (Q - Qalt) over the element, where Qalt is the truncated version of Q
		Here we don't bother using quadrature, we do a simple sum
	*/
	UClipped = sf.Clipper.Mul(U, UClipped)
	var mNum, mDenom float64
	for i := 0; i < Np; i++ {
		mass := MD[i+i*Np]
		t1 := UD[i] - UClippedD[i]
		mNum += mass * (t1 * t1)
		mDenom += mass * (UD[i] * UD[i])
	}
	m = mNum / mDenom
	return
}

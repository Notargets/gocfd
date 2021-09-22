package Euler2D

import (
	"math"

	"github.com/notargets/gocfd/DG2D"
	"github.com/notargets/gocfd/utils"
)

type ModeAliasShockFinder struct {
	Element *DG2D.LagrangeElement2D
	Clipper utils.Matrix // Matrix used to clip the topmost mode from the solution polynomial, used in shockfinder
	Np      int
	q, qalt utils.Matrix // scratch storage for evaluating the moment
}

func NewAliasShockFinder(dfr *DG2D.LagrangeElement2D) (sf *ModeAliasShockFinder) {
	var (
		Np = dfr.Np
	)
	sf = &ModeAliasShockFinder{
		Element: dfr,
		Np:      Np,
		q:       utils.NewMatrix(Np, 1),
		qalt:    utils.NewMatrix(Np, 1),
	}
	data := make([]float64, Np)
	for i := 0; i < Np; i++ {
		if i != Np-1 {
			data[i] = 1.
		} else {
			data[i] = 0.
		}
	}
	diag := utils.NewDiagMatrix(Np, data)
	/*
		The "Clipper" matrix drops the last mode from the polynomial and forms an alternative field of values at the node
		points based on a polynomial with one less term. In other words, if we have a polynomial of degree "p", expressed
		as values at Np node points, multiplying the Node point values vector by Clipper produces an alternative version
		of the node values based on truncating the last polynomial mode.
	*/
	sf.Clipper = dfr.V.Mul(diag).Mul(dfr.Vinv)
	return
}

func (sf *ModeAliasShockFinder) ShockIndicator(q []float64) (sigma float64) {
	/*
		Original method by Perrson, constants chosen to match Zhiqiang, et. al.
		Zhiqiang uses a threshold of sigma<0.99 to indicate "troubled cell"
	*/
	var (
		Se          = math.Log10(sf.moment(q))
		k           = float64(sf.Element.N)
		kappa       = 4.
		C0          = 3.
		S0          = -C0 * math.Log(k)
		left, right = S0 - kappa, S0 + kappa
		ookappa     = 1. / kappa
	)
	switch {
	case Se < left:
		sigma = 1.
	case Se >= left && Se < right:
		sigma = 0.5 * (1. - math.Sin(0.5*math.Pi*ookappa*(Se-S0)))
	case Se >= right:
		sigma = 0.
	}
	return
}

func (sf *ModeAliasShockFinder) moment(q []float64) (m float64) {
	var (
		qd, qaltd = sf.q.DataP, sf.qalt.DataP
	)
	/*
		Evaluate the L2 moment of (q - qalt) over the element, where qalt is the truncated version of q
		Here we don't bother using quadrature, we do a simple sum
	*/
	copy(sf.q.DataP, q)
	sf.qalt = sf.Clipper.Mul(sf.q, sf.qalt)
	for i := 0; i < sf.Np; i++ {
		t1 := qd[i] - qaltd[i]
		m += t1 * t1 / (qd[i] * qd[i])
	}
	return
}

package Euler2D

import (
	"fmt"
	"math"

	"github.com/notargets/gocfd/DG2D"
	"github.com/notargets/gocfd/utils"
)

type BarthJespersonLimiter struct {
	Element     *DG2D.LagrangeElement2D
	Tris        *DG2D.Triangulation
	ShockFinder *ModeAliasShockFinder
	UElement    []float64 // Scratch area for assembly and testing of solution values
}

func NewBarthJespersonLimiter(dfr DG2D.DFR2D, sf *ModeAliasShockFinder) (bjl *BarthJespersonLimiter) {
	bjl = &BarthJespersonLimiter{
		Element:     dfr.SolutionElement,
		Tris:        dfr.Tris,
		ShockFinder: sf,
		UElement:    make([]float64, dfr.SolutionElement.Np),
	}
	return
}

func (bjl *BarthJespersonLimiter) LimitSolution(Q [4]utils.Matrix) {
	var (
		Np, Kmax         = Q[0].Dims()
		Uave, Umin, Umax float64
	)
	getElAvg := func(f utils.Matrix, k int) (ave float64) {
		for i := 0; i < Np; i++ {
			ind := k + Kmax*i
			ave += f.DataP[ind]
		}
		ave /= float64(Np)
		return
	}
	for k := 0; k < Kmax; k++ {
		for i := 0; i < Np; i++ {
			ind := k + Kmax*i
			bjl.UElement[i] = Q[3].DataP[ind] // Use Energy as the indicator basis
		}
		if bjl.ShockFinder.ElementHasShock(bjl.UElement) { // Element has a shock
			// Apply limiting procedure
			// Get average solution value for element and neighbors
			Uave = getElAvg(Q[3], k)
			Umin, Umax = Uave, Uave
			// Loop over connected tris to get Umin, Umax
			for ii := 0; ii < 3; ii++ {
				kk := bjl.Tris.EtoE[k][ii]
				if kk != -1 {
					U := getElAvg(Q[3], kk)
					Umax = math.Max(U, Umax)
					Umin = math.Min(U, Umin)
				}
			}
			// Obtain average gradient of this cell
		}
	}
}

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

func (sf *ModeAliasShockFinder) ElementHasShock(q []float64) (i bool) {
	// Zhiqiang uses a threshold of sigma<0.99 to indicate "troubled cell"
	if sf.ShockIndicator(q) < 0.99 {
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
		k           = float64(sf.Element.N)
		kappa       = 4.
		C0          = 3.
		S0          = -C0 * math.Log10(k)
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
	if len(q) != sf.Np {
		err := fmt.Errorf("incorrect dimension of solution vector, should be %d is %d",
			sf.Np, len(q))
		panic(err)
	}
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

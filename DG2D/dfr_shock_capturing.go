package DG2D

import (
	"math"

	"github.com/notargets/gocfd/utils"
)

type ModeAliasShockFinder struct {
	Element             *LagrangeElement2D
	Clipper             utils.Matrix // Matrix used to clip the topmost mode from the solution polynomial, used in shockfinder
	P, Mass, D          utils.Matrix // Element matrices
	Np                  int
	Q, Qalt             utils.Matrix // scratch storage for evaluating the moment
	Kappa               float64
	ShockCells          *utils.DynBuffer[int]
	ShockSigmaThreshold float64
}

func (sf *ModeAliasShockFinder) UpdateShockedCells(Rho utils.Matrix) {
	var (
		_, KMax = Rho.Dims()
	)
	if sf.ShockCells == nil {
		sf.ShockCells = utils.NewDynBuffer[int](KMax)
	}
	sf.ShockCells.Reset()
	for k := 0; k < KMax; k++ {
		scratch := sf.Qalt.DataP
		for i := 0; i < sf.Np; i++ {
			ind := k + i*KMax
			scratch[i] = Rho.DataP[ind]
		}
		if sf.ElementHasShock(scratch) {
			sf.ShockCells.Add(k)
		}
	}
	return
}

func (dfr *DFR2D) NewAliasShockFinder(Kappa float64) (sf *ModeAliasShockFinder) {
	var (
		element = dfr.SolutionElement
		Np      = element.Np
		N       = element.N
	)
	sf = &ModeAliasShockFinder{
		Element:             element,
		Np:                  Np,
		Q:                   utils.NewMatrix(Np, 1),
		Qalt:                utils.NewMatrix(Np, 1),
		Kappa:               Kappa,
		ShockSigmaThreshold: 0.0075,
		Mass:                element.MassMatrix,
	}
	/*
		The "Clipper" matrix drops the last mode from the polynomial and forms an alternative field of values at the node
		points based on a polynomial with one less term. In other words, if we have a polynomial of degree "p", expressed
		as values at Np node points, multiplying the Node point values vector by Clipper produces an alternative version
		of the node values based on truncating the last polynomial mode.
	*/
	// Compute element matrices
	sf.Clipper =
		element.JB2D.V.Mul(dfr.CutoffFilter2D(N, N, 0)).Mul(element.JB2D.Vinv)

	// sf.D = utils.NewMatrix(Np, Np).AddScalar(1.).Subtract(sf.Clipper)
	sf.D = utils.NewDiagMatrix(Np, nil, 1.).Subtract(sf.Clipper)
	sf.P = sf.Mass.Mul(sf.D)
	return
}

func (sf *ModeAliasShockFinder) UpdateSeMoment(TestVar utils.Matrix,
	LScratch [3]utils.Matrix, Se utils.Vector) {
	// X and Y are scratch matrices of size Np,KMax
	// Se should be a vector of size KMax
	var (
		Np, KMax = TestVar.Dims()
		X, Y, DQ = LScratch[0], LScratch[1], LScratch[2]
	)
	X = sf.P.Mul(TestVar, X)
	Y = sf.Mass.Mul(TestVar, Y)
	DQ = sf.D.Mul(TestVar, DQ)
	for k := 0; k < KMax; k++ {
		var num, den float64
		for i := 0; i < Np; i++ {
			// d_i = (U - Ualt)[i,k]
			// but Ualt = C*Q so U - Ualt = D*Q, and X = M * D * Q already
			di := DQ.At(i, k)
			num += di * X.At(i, k)
			den += TestVar.At(i, k) * Y.At(i, k)
		}
		Se.Set(k, math.Log10(num/den))
	}
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
	// if sf.ShockIndicator(q) > 0.01 {
	if sf.ShockIndicator(q) > sf.ShockSigmaThreshold {
		i = true
	}
	return
}

func (sf *ModeAliasShockFinder) GetShockIndicator(Q utils.Matrix, k int) (sigma float64) {
	var (
		Np, Kmax = Q.Dims()
		Dens     = Q.DataP
		scratch  = sf.Qalt.DataP
	)
	for i := 0; i < Np; i++ {
		ind := k + i*Kmax
		scratch[i] = Dens[ind]
	}
	return sf.ShockIndicator(scratch)
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

		Note: The mass matrix of the basis is diagonal, as it is orthogonal
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

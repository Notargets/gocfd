package DG2D

import (
	"fmt"
	"math"

	"github.com/notargets/gocfd/utils"
)

type ModeAliasShockFinder struct {
	Element             *LagrangeElement2D
	Clipper             utils.Matrix // Matrix used to clip the topmost mode from the solution polynomial, used in shockfinder
	ModeFilter          []float64
	P, MassMatrix, D    utils.Matrix // Element matrices
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

func (dfr *DFR2D) ModalFilter(alpha float64, s int) (mf []float64) {
	// This should be used in conjunction with the limiter,
	// then applied to the modes directly.
	var (
		N  = dfr.SolutionElement.N
		Np = dfr.SolutionElement.Np
	)
	mf = make([]float64, Np)

	if N <= 1 {
		for j := range mf {
			mf[j] = 1.
		}
	} else {
		for j := range mf {
			p := dfr.SolutionElement.JB2D.OrderAtJ[j]
			if p == 0 {
				mf[j] = 1.0
			} else {
				normalized := float64(p) / float64(N)
				// Interpolate between 1 (no damping) and baseDamp (max damping)
				mf[j] = math.Exp(-alpha * math.Pow(normalized, float64(s)))
			}
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
		MassMatrix:          element.MassMatrix,
	}
	// Compute element matrices
	/*
		The "Clipper" matrix drops the last mode from the polynomial and forms an alternative field of values at the node
		points based on a polynomial with one less term. In other words, if we have a polynomial of degree "p", expressed
		as values at Np node points, multiplying the Node point values vector by Clipper produces an alternative version
		of the node values based on truncating the last polynomial mode.
	*/
	// sf.Clipper =
	// 	element.JB2D.V.Mul(dfr.CutoffFilter2D(N, 0)).Mul(element.JB2D.Vinv)
	sf.Clipper =
		element.JB2D.V.Mul(dfr.CutoffFilter2D(N, 0)).Mul(element.JB2D.Vinv)
	// Implement a cutoff filter to suppress ringing, won't alter modes < 3
	alpha, s := RecommendedFilterParameters(dfr.N)
	sf.ModeFilter = dfr.ModalFilter(alpha, s)
	// fmt.Println("Exponential Filter: ", sf.ModeFilter)
	// os.Exit(1)

	// sf.D = utils.NewMatrix(Np, Np).AddScalar(1.).Subtract(sf.Clipper)
	sf.D = utils.NewDiagMatrix(Np, nil, 1.).Subtract(sf.Clipper)
	sf.P = sf.MassMatrix.Mul(sf.D)
	return
}

func RecommendedFilterParameters(P int) (alpha float64, s int) {
	switch {
	case P <= 2:
		alpha = 4.0
		s = 4
	case P == 3:
		alpha = 6.0
		s = 4
	case P == 4:
		alpha = 8.0
		s = 5
	case P == 5:
		alpha = 10.0
		s = 6
		// s = 48
	case P == 6:
		alpha = 14.0
		s = 6
	case P == 7:
		alpha = 18.0
		s = 7
	case P == 8:
		alpha = 24.0
		s = 8
	case P >= 9:
		alpha = 30.0 + 2.0*float64(P-9)
		s = 8 + (P-9)/2
	default:
		panic(fmt.Errorf("invalid polynomial order P = %d", P))
	}
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
	Y = sf.MassMatrix.Mul(TestVar, Y)
	DQ = sf.D.Mul(TestVar, DQ)
	for k := 0; k < KMax; k++ {
		var num, den float64
		for i := 0; i < Np; i++ {
			ind := k + i*KMax
			// d_i = (U - Ualt)[i,k]
			// but Ualt = C*Q so U - Ualt = D*Q, and X = M * D * Q already
			di := DQ.DataP[ind]
			num += di * X.DataP[ind]
			den += TestVar.DataP[ind] * Y.DataP[ind]
		}
		Se.Set(k, math.Log10(num/den))
	}
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

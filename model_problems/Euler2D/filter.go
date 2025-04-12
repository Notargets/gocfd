package Euler2D

import (
	"fmt"
	"math"
	"strings"

	"github.com/notargets/gocfd/DG2D"
	"github.com/notargets/gocfd/utils"
)

type SolutionLimiter struct {
	limiterType          LimiterType
	Element              *DG2D.LagrangeElement2D
	Tris                 *DG2D.Triangulation
	Partitions           *utils.PartitionMap
	ShockFinder          []*DG2D.ModeAliasShockFinder // Sharded
	UElement, dUdr, dUds []utils.Matrix               // Sharded scratch areas for assembly and testing of solution values
	FS                   *FreeStream
}

type LimiterType uint8

const (
	None LimiterType = iota
	BarthJespersonT
	PerssonC0T
)

var (
	LimiterNames = map[string]LimiterType{
		"barthjesperson":  BarthJespersonT,
		"barth jesperson": BarthJespersonT,
		"perssonc0":       PerssonC0T,
		"persson c0":      PerssonC0T,
	}
	LimiterNamesRev = map[LimiterType]string{
		BarthJespersonT: "Barth Jesperson",
		PerssonC0T:      "Persson, C0 viscosity",
	}
)

func (lt LimiterType) Print() (txt string) {
	if val, ok := LimiterNamesRev[lt]; !ok {
		txt = "None"
	} else {
		txt = val
	}
	return
}

func NewLimiterType(label string) (lt LimiterType) {
	var (
		ok  bool
		err error
	)
	if len(label) == 0 {
		return None
	}
	label = strings.ToLower(strings.TrimSpace(label))
	if lt, ok = LimiterNames[label]; !ok {
		err = fmt.Errorf("unable to use limiter named [%s]", label)
		panic(err)
	}
	return
}

func NewSolutionLimiter(t LimiterType, kappa float64, dfr *DG2D.DFR2D, pm *utils.PartitionMap, fs *FreeStream) (bjl *SolutionLimiter) {
	var (
		Np       = dfr.SolutionElement.Np
		Nthreads = pm.ParallelDegree
	)
	bjl = &SolutionLimiter{
		limiterType: t,
		Element:     dfr.SolutionElement,
		Tris:        dfr.Tris,
		ShockFinder: make([]*DG2D.ModeAliasShockFinder, Nthreads),
		FS:          fs,
		Partitions:  pm,
		// Sharded working matrices
		UElement: make([]utils.Matrix, Nthreads),
		dUdr:     make([]utils.Matrix, Nthreads),
		dUds:     make([]utils.Matrix, Nthreads),
	}
	switch bjl.limiterType {
	case BarthJespersonT:
		for np := 0; np < Nthreads; np++ {
			bjl.ShockFinder[np] = dfr.NewAliasShockFinder(kappa)
			bjl.UElement[np] = utils.NewMatrix(Np, 1)
			bjl.dUdr[np] = utils.NewMatrix(Np, 1)
			bjl.dUds[np] = utils.NewMatrix(Np, 1)
		}
	}
	return
}

func (bjl *SolutionLimiter) LimitSolution(myThread int, Qall, Residual [][4]utils.Matrix, FilterScratch []utils.Matrix) (points int) {
	var (
		Q        = Qall[myThread]
		Np, Kmax = Q[0].Dims()
		UE       = bjl.UElement[myThread]
	)
	switch bjl.limiterType {
	case BarthJespersonT:
		for k := 0; k < Kmax; k++ {
			for i := 0; i < Np; i++ {
				ind := k + Kmax*i
				UE.DataP[i] = Q[0].DataP[ind]
			}
			if bjl.ShockFinder[myThread].ElementHasShock(UE.DataP) { // Element has a shock
				bjl.limitScalarFieldBarthJesperson(k, myThread, Qall)
				for n := 0; n < 4; n++ {
					for i := 0; i < Np; i++ {
						ind := k + Kmax*i
						Residual[myThread][n].DataP[ind] = 0.
					}
				}
				points++
			}
		}
	}
	return
}

func (bjl *SolutionLimiter) limitScalarFieldBarthJesperson(k, myThread int, Qall [][4]utils.Matrix) {
	var (
		el         = bjl.Element
		Np, Kmax   = Qall[myThread][0].Dims()
		Dr, Ds     = el.Dr, el.Ds
		MMD        = el.MassMatrix.DataP
		UE         = bjl.UElement[myThread]
		dUdr, dUds = bjl.dUdr[myThread], bjl.dUds[myThread]
		min, max   = math.Min, math.Max
	)
	getElAvg := func(f utils.Matrix, kkk, kMx int) (ave float64) {
		var massTotal float64
		for i := 0; i < Np; i++ {
			ind := kkk + kMx*i
			mass := MMD[i+i*Np]
			massTotal += mass
			ave += mass * f.DataP[ind]
		}
		ave /= massTotal
		return
	}
	psiCalc := func(corner int, Umin, Umax, Uave, dUdrAve, dUdsAve float64) (psi float64) {
		var (
			del2 float64
		)
		/*
			For each corner of the unit triangle, the vector from center to corner is:
				ri[0] = [ -2/3, -2/3 ]
				ri[1] = [ 4/3, -2/3 ]
				ri[2] = [ -2/3, 4/3 ]
		*/
		switch corner {
		case 0:
			del2 = -2. / 3. * (dUdrAve + dUdsAve)
		case 1:
			del2 = (4./3.)*dUdrAve - (2./3.)*dUdsAve
		case 2:
			del2 = -(2./3.)*dUdrAve + (4./3.)*dUdsAve
		}
		oodel2 := 1. / del2
		// Calculate limiter function Psi
		switch {
		case del2 > 0:
			psi = min(1, oodel2*(Umax-Uave))
		case del2 == 0:
			psi = 1
		case del2 < 0:
			psi = min(1, oodel2*(Umin-Uave))
		}
		return
	}
	for n := 0; n < 4; n++ {
		var (
			U                = Qall[myThread][n]
			Uave, Umin, Umax float64
		)
		// Apply limiting procedure
		// Get average and min/max solution value for element and neighbors
		Uave = getElAvg(U, k, Kmax)
		Umin, Umax = Uave, Uave
		// Loop over connected tris to get Umin, Umax
		for ii := 0; ii < 3; ii++ {
			kk := bjl.Tris.EtoE[k][ii]
			if kk != -1 {
				remoteK, remoteKmax, rThread := bjl.Partitions.GetLocalK(kk)
				UU := getElAvg(Qall[rThread][n], remoteK, remoteKmax)
				Umax = max(UU, Umax)
				Umin = min(UU, Umin)
			}
		}
		for i := 0; i < Np; i++ {
			ind := k + Kmax*i
			UE.DataP[i] = U.DataP[ind]
		}
		// Obtain average gradient of this cell
		getAvgDeriv := func(deriv utils.Matrix) (derivAve float64) {
			var (
				massTotal float64
				derivD    = deriv.DataP
			)
			for i := 0; i < Np; i++ {
				mass := MMD[i+i*Np]
				massTotal += mass
				derivAve += mass * derivD[i]
			}
			derivAve /= massTotal
			return
		}
		dUdrAve, dUdsAve := getAvgDeriv(Dr.Mul(UE, dUdr)), getAvgDeriv(Ds.Mul(UE, dUds))
		// Form psi as the minimum of all three corners
		var psi float64
		psi = 10000.
		for nn := 0; nn < 3; nn++ {
			psi = min(psi, psiCalc(nn, Umin, Umax, Uave, dUdrAve, dUdsAve))
		}
		// Limit the solution using psi and the average gradient
		for i := 0; i < Np; i++ {
			// Vector from node points to center of element
			dR, dS := bjl.Element.R.DataP[i]-(-1./3), bjl.Element.S.DataP[i]-(-1./3.)
			ind := k + Kmax*i
			U.DataP[ind] = Uave + psi*(dR*dUdrAve+dS*dUdsAve)
		}
	}
}

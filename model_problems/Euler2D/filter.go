package Euler2D

import (
	"fmt"
	"math"
	"strings"

	"github.com/notargets/gocfd/DG2D"
	"github.com/notargets/gocfd/utils"
)

type LimiterType uint8

const (
	None LimiterType = iota
	PerssonC0T
)

var (
	LimiterNames = map[string]LimiterType{
		"perssonc0":  PerssonC0T,
		"persson c0": PerssonC0T,
	}
	LimiterNamesRev = map[LimiterType]string{
		PerssonC0T: "Persson, C0 viscosity",
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

func LimitSolution(Q [4]utils.Matrix, QMean [4]utils.Vector,
	Sigma utils.Vector, sf *DG2D.ModeAliasShockFinder) (
	points int) {
	var (
		Np, Kmax = Q[0].Dims()
		Beta     = 10.
	)
	for k := 0; k < Kmax; k++ {
		sigma := Sigma.AtVec(k)
		// if sigma > sf.ShockSigmaThreshold { // Element has a shock
		alpha := 1. - math.Exp(-Beta*sigma)
		// alpha := math.Pow(sigma, 2.)
		// alpha := sigma
		// alpha := math.Pow(sigma, 1./3.)
		// fmt.Printf("sigma, alpha[%d] = %.1f, %.1f\n", k, sigma, alpha)
		for n := 0; n < 4; n++ {
			for i := 0; i < Np; i++ {
				Q[n].Set(i, k,
					(1.-alpha)*Q[n].At(i, k)+alpha*QMean[n].AtVec(k))
			}
		}
		// }
	}
	return
}

package Euler2D

import (
	"fmt"
	"math"
	"strings"

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

func LimitSolution(myThread int, Q [4]utils.Matrix, QMean [4]utils.Vector,
	sd *ScalarDissipation) (
	points int) {
	var (
		Np, Kmax    = Q[0].Dims()
		SigmaScalar = sd.SigmaScalar[myThread]
		Sigma       = sd.Sigma[myThread]
	)
	switch sd.Continuity {
	case No:
		for k := 0; k < Kmax; k++ {
			sigma := SigmaScalar.AtVec(k)
			alpha := math.Pow(sigma, 1./2.)
			for n := 0; n < 4; n++ {
				for i := 0; i < Np; i++ {
					Q[n].Set(i, k,
						(1.-alpha)*Q[n].At(i, k)+alpha*QMean[n].AtVec(k))
				}
			}
		}
	case C0:
		for n := 0; n < 4; n++ {
			for k := 0; k < Kmax; k++ {
				for i := 0; i < Np; i++ {
					sigma := Sigma.At(i, k)
					alpha := math.Pow(sigma, 1./2.)
					Q[n].Set(i, k,
						(1.-alpha)*Q[n].At(i, k)+alpha*QMean[n].AtVec(k))
				}
			}
		}
	}
	return
}

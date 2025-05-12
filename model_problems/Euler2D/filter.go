package Euler2D

import (
	"fmt"
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

func LimitSolution(Q [4]utils.Matrix,
	Sigma utils.Vector, sf *DG2D.ModeAliasShockFinder) (points int) {
	var (
		Np, Kmax = Q[0].Dims()
	)
	_ = Np
	for k := 0; k < Kmax; k++ {
		if Sigma.AtVec(k) > sf.ShockSigmaThreshold { // Element has a shock
			// TODO: Implement sigma limiter
		}
	}
	return
}

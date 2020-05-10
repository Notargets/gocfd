package main

import (
	"flag"
	"fmt"
	"math"
	"time"

	"github.com/notargets/gocfd/model_problems"
)

var (
	K         = 0 // Number of elements
	N         = 0 // Polynomial degree
	Delay     = time.Duration(0)
	ModelRun  = Euler1D
	CFL       = 0.0
	FinalTime = 100000.
	XMax      = 0.0
)

type ModelType uint8

const (
	Advect1D ModelType = iota
	Maxwell1D
	Euler1D
	AdvectDFR
	MaxwellDFR
)

type Model interface {
	Run(graph bool, graphDelay ...time.Duration)
}

func main() {
	ModelRunptr := flag.Int("model", int(ModelRun), "model to run: 0 = Advect1D, 1 = Maxwell1D, 2 = Euler1D")
	Kptr := flag.Int("K", K, "Number of elements in model")
	Nptr := flag.Int("N", N, "polynomial degree")
	Delayptr := flag.Int("delay", 0, "milliseconds of delay for plotting")
	Graphptr := flag.Bool("graph", false, "display a graph while computing solution")
	CFLptr := flag.Float64("CFL", CFL, "CFL - increase for speedup, decrease for stability")
	FTptr := flag.Float64("FinalTime", FinalTime, "FinalTime - the target end time for the sim")
	XMaxptr := flag.Float64("XMax", XMax, "Maximum X coordinate (for Euler) - make sure to increase K with XMax")
	flag.Parse()
	ModelRun = ModelType(*ModelRunptr)
	Delay = time.Duration(*Delayptr)
	CFL, XMax, N, K = Defaults(ModelRun)

	K = int(getParam(float64(K), Kptr))
	N = int(getParam(float64(N), Nptr))
	CFL = LimitCFL(ModelRun, getParam(CFL, CFLptr))
	FinalTime = getParam(FinalTime, FTptr)
	XMax = getParam(XMax, XMaxptr)

	var C Model
	switch ModelRun {
	case Advect1D:
		C = model_problems.NewAdvection1D(2*math.Pi, CFL, FinalTime, N, K)
	case Maxwell1D:
		C = model_problems.NewMaxwell1D(CFL, FinalTime, N, K)
	case AdvectDFR:
		C = model_problems.NewAdvectionDFR(2*math.Pi, CFL, FinalTime, XMax, N, K)
	case MaxwellDFR:
		C = model_problems.NewMaxwellDFR(CFL, FinalTime, N, K)
	case Euler1D:
		fallthrough
	default:
		C = model_problems.NewEuler1D(CFL, FinalTime, XMax, N, K)
	}
	C.Run(*Graphptr, Delay*time.Millisecond)
}

func LimitCFL(ModelRun ModelType, CFL float64) (CFLNew float64) {
	var (
		CFLMax float64
	)
	switch ModelRun {
	case Advect1D:
		CFLMax = 1
	case Maxwell1D:
		CFLMax = 1
	case Euler1D:
		CFLMax = 3
	case AdvectDFR:
		CFLMax = 3
	case MaxwellDFR:
		CFLMax = 1
	}
	if CFL > CFLMax {
		fmt.Printf("Input CFL is higher than max CFL for this method\nReplacing with Max CFL: %8.2f\n", CFLMax)
		return CFLMax
	}
	return CFL
}

func Defaults(ModelRun ModelType) (CFL, XMax float64, N, K int) {
	switch ModelRun {
	case Advect1D:
		K = 10
		N = 3
		CFL = 1
		XMax = 2 * math.Pi
	case Maxwell1D:
		K = 100
		N = 4
		CFL = 1
		XMax = 1
	case Euler1D:
		K = 500
		N = 4
		CFL = 3
		XMax = 1
	case AdvectDFR:
		K = 50
		N = 4
		CFL = 3
		XMax = 2 * math.Pi
	case MaxwellDFR:
		K = 100
		N = 4
		CFL = 1
		XMax = 1
	}
	return
}

func getParam(def float64, valP interface{}) float64 {
	switch val := valP.(type) {
	case *int:
		if *val != 0 {
			return float64(*val)
		}
	case *float64:
		if *val != 0 {
			return *val
		}
	}
	return def
}

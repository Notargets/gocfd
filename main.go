package main

import (
	"flag"
	"fmt"
	"math"
	"time"

	"github.com/notargets/gocfd/model_problems/Advection1D"
	Maxwell1D2 "github.com/notargets/gocfd/model_problems/Maxwell1D"

	Euler1D2 "github.com/notargets/gocfd/model_problems/Euler1D"
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
	EulerDFR_Roe
	EulerDFR_LF
)

var (
	max_CFL  = []float64{1, 1, 3, 3, 1, 2.2, 3}
	def_K    = []int{10, 100, 500, 50, 500, 500, 500}
	def_N    = []int{3, 4, 4, 4, 3, 4, 4}
	def_CFL  = []float64{1, 1, 3, 3, 0.75, 2.2, 3}
	def_XMAX = []float64{2 * math.Pi, 1, 1, 2 * math.Pi, 1, 1, 1}
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
		C = Advection1D.NewAdvection1D(2*math.Pi, CFL, FinalTime, N, K)
	case Maxwell1D:
		C = Maxwell1D2.NewMaxwellDFR(CFL, FinalTime, N, K, Maxwell1D2.GK)
	case AdvectDFR:
		C = Advection1D.NewAdvectionDFR(2*math.Pi, CFL, FinalTime, XMax, N, K)
	case MaxwellDFR:
		C = Maxwell1D2.NewMaxwellDFR(CFL, FinalTime, N, K, Maxwell1D2.DFR)
	case EulerDFR_Roe:
		C = Euler1D2.NewEulerDFR(CFL, FinalTime, XMax, N, K, Euler1D2.Euler_DFR_Roe)
	case EulerDFR_LF:
		C = Euler1D2.NewEulerDFR(CFL, FinalTime, XMax, N, K, Euler1D2.Euler_DFR_LF)
	case Euler1D:
		fallthrough
	default:
		//C = model_problems.NewEuler1D(CFL, FinalTime, XMax, N, K)
		C = Euler1D2.NewEulerDFR(CFL, FinalTime, XMax, N, K, Euler1D2.Galerkin_LF)
	}
	C.Run(*Graphptr, Delay*time.Millisecond)
}

func LimitCFL(model ModelType, CFL float64) (CFLNew float64) {
	var (
		CFLMax float64
	)
	CFLMax = max_CFL[model]
	if CFL > CFLMax {
		fmt.Printf("Input CFL is higher than max CFL for this method\nReplacing with Max CFL: %8.2f\n", CFLMax)
		return CFLMax
	}
	return CFL
}

func Defaults(model ModelType) (CFL, XMax float64, N, K int) {
	return def_CFL[model], def_XMAX[model], def_N[model], def_K[model]
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

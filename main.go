package main

import (
	"flag"
	"math"
	"time"

	"github.com/notargets/gocfd/model_problems"
)

var (
	K         = 60 // Number of elements
	N         = 8  // Polynomial degree
	Delay     = time.Duration(0)
	ModelRun  = Maxwell1D
	CFL       = 1.0
	FinalTime = 100000.
)

type ModelType uint8

const (
	Advect1D ModelType = iota
	Maxwell1D
	Euler1D
)

type Model interface {
	Run(graph bool, graphDelay ...time.Duration)
}

func main() {
	Kptr := flag.Int("K", K, "Number of elements in model")
	Nptr := flag.Int("N", N, "polynomial degree")
	Delayptr := flag.Int("delay", 0, "milliseconds of delay for plotting")
	Graphptr := flag.Bool("graph", false, "display a graph while computing solution")
	ModelRunptr := flag.Int("model", int(ModelRun), "model to run: 0 = Advect1D, 1 = Maxwell1D, 2 = Euler1D")
	CFLptr := flag.Float64("CFL", CFL, "CFL - increase for speedup, decrease for stability")
	FTptr := flag.Float64("FinalTime", FinalTime, "FinalTime - the target end time for the sim")
	flag.Parse()
	K = *Kptr
	N = *Nptr
	Delay = time.Duration(*Delayptr)
	ModelRun = ModelType(*ModelRunptr)
	CFL = *CFLptr
	FinalTime = *FTptr

	var C Model
	switch ModelRun {
	case Advect1D:
		C = model_problems.NewAdvection1D(2*math.Pi, CFL, FinalTime, N, K)
	case Euler1D:
		C = model_problems.NewEuler1D(CFL, FinalTime, N, K)
	case Maxwell1D:
		fallthrough
	default:
		C = model_problems.NewMaxwell1D(CFL, FinalTime, N, K)
	}
	C.Run(*Graphptr, Delay*time.Millisecond)
}

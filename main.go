package main

import (
	"flag"
	"math"
	"time"

	"github.com/notargets/gocfd/model_problems"
)

var (
	K        = 60 // Number of elements
	N        = 8  // Polynomial degree
	Delay    = time.Duration(0)
	ModelRun = Maxwell1D
)

type ModelType uint8

const (
	Advect1D ModelType = iota
	Maxwell1D
)

type Model interface {
	Run(graph bool, graphDelay ...time.Duration)
}

func main() {
	Kptr := flag.Int("K", K, "Number of elements in model")
	Nptr := flag.Int("N", N, "polynomial degree")
	Delayptr := flag.Int("delay", 200, "milliseconds of delay for plotting")
	Graphptr := flag.Bool("graph", false, "display a graph while computing solution")
	ModelRunptr := flag.Int("model", int(ModelRun), "model to run: 0 = Advect1D, 1 = Maxwell1D")
	flag.Parse()
	K = *Kptr
	N = *Nptr
	Delay = time.Duration(*Delayptr)
	ModelRun = ModelType(*ModelRunptr)

	var C Model
	switch ModelRun {
	case Advect1D:
		C = model_problems.NewAdvection1D(2*math.Pi, 0.75, 100000., N, K)
	case Maxwell1D:
		fallthrough
	default:
		C = model_problems.NewMaxwell1D(1.0, 100000., N, K)
	}
	C.Run(*Graphptr, Delay*time.Millisecond)
}

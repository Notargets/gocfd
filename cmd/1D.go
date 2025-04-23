/*
Copyright © 2020 NAME HERE <EMAIL ADDRESS>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

	http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/
package cmd

import (
	"fmt"
	"math"
	"time"

	"github.com/notargets/gocfd/model_problems/Advection1D"
	"github.com/notargets/gocfd/model_problems/Euler1D"
	"github.com/notargets/gocfd/model_problems/Maxwell1D"
	"github.com/spf13/cobra"
)

// OneDCmd represents the 1D command
var OneDCmd = &cobra.Command{
	Use:   "1D",
	Short: "One Dimensional Model Problem Solutions",
	Long: `
Executes the Nodal Discontinuous Galerkin solver for a variety of model problems,

gocfd 1D `,
	Run: func(cmd *cobra.Command, args []string) {
		m1d := &Model1D{}
		fmt.Println("1D called")
		mr, _ := cmd.Flags().GetInt("model")
		m1d.ModelRun = ModelType1D(mr)
		Casep, _ := cmd.Flags().GetInt("case")
		m1d.Case = Euler1D.CaseType(Casep)
		m1d.XMax, _ = cmd.Flags().GetFloat64("xMax")
		m1d.FinalTime, _ = cmd.Flags().GetFloat64("finalTime")
		m1d.CFL, _ = cmd.Flags().GetFloat64("CFL")
		m1d.N, _ = cmd.Flags().GetInt("n")
		m1d.K, _ = cmd.Flags().GetInt("k")
		m1d.CFL = LimitCFL(m1d.ModelRun, m1d.CFL)
		Run1D(m1d)
	},
}

func init() {
	rootCmd.AddCommand(OneDCmd)
	var (
		K         = 0 // Number of elements
		N         = 0 // Polynomial degree
		ModelRun  = M_1DEuler
		CFL       = 0.0
		FinalTime = 100000.
		XMax      = 0.0
		CaseInt   int
	)
	CFL, XMax, N, K, CaseInt = Defaults(ModelRun)
	OneDCmd.Flags().IntP("model", "m", int(ModelRun), "model to run: 0 = Advect1D, 1 = Maxwell1D, 2 = Euler1D")
	OneDCmd.Flags().IntP("k", "k", K, "Number of elements in model")
	OneDCmd.Flags().IntP("n", "n", N, "polynomial degree")
	OneDCmd.Flags().IntP("delay", "d", 0, "milliseconds of delay for plotting")
	OneDCmd.Flags().IntP("case", "c", int(CaseInt), "Case to run, for Euler: 0 = SOD Shock Tube, 1 = Density Wave")
	OneDCmd.Flags().Float64("CFL", CFL, "CFL - increase for speedup, decrease for stability")
	OneDCmd.Flags().Float64("finalTime", FinalTime, "FinalTime - the target end time for the sim")
	OneDCmd.Flags().Float64("xMax", XMax, "Maximum R coordinate (for Euler) - make sure to increase K with XMax")
}

type Model1D struct {
	K, N                 int // Number of elements, Polynomial Degree
	Delay                time.Duration
	ModelRun             ModelType1D
	CFL, FinalTime, XMax float64
	Case                 Euler1D.CaseType
	Graph                bool
}

type ModelType1D uint8

const (
	M_1DAdvect ModelType1D = iota
	M_1DMaxwell
	M_1DEuler
	M_1DAdvectDFR
	M_1DMaxwellDFR
	M_1DEulerDFR_Roe
	M_1DEulerDFR_LF
	M_1DEulerDFR_Ave
)

var (
	max_CFL  = []float64{1, 1, 3, 3, 1, 2.5, 3, 3}
	def_K    = []int{10, 100, 500, 50, 500, 500, 500, 40}
	def_N    = []int{3, 4, 4, 4, 3, 4, 4, 3}
	def_CFL  = []float64{1, 1, 3, 3, 0.75, 2.5, 3, 0.5}
	def_XMAX = []float64{2 * math.Pi, 1, 1, 2 * math.Pi, 1, 1, 1, 1}
	def_CASE = make([]int, 8)
)

type Model interface {
	Run(graph bool, graphDelay ...time.Duration)
}

func Run1D(m1d *Model1D) {
	var C Model
	switch m1d.ModelRun {
	case M_1DAdvect:
		C = Advection1D.NewAdvection(2*math.Pi, m1d.CFL, m1d.FinalTime, m1d.XMax, m1d.N, m1d.K, Advection1D.GK)
	case M_1DMaxwell:
		C = Maxwell1D.NewMaxwell(m1d.CFL, m1d.FinalTime, m1d.N, m1d.K, Maxwell1D.GK)
	case M_1DAdvectDFR:
		C = Advection1D.NewAdvection(2*math.Pi, m1d.CFL, m1d.FinalTime, m1d.XMax, m1d.N, m1d.K, Advection1D.DFR)
	case M_1DMaxwellDFR:
		C = Maxwell1D.NewMaxwell(m1d.CFL, m1d.FinalTime, m1d.N, m1d.K, Maxwell1D.DFR)
	case M_1DEulerDFR_Roe:
		C = Euler1D.NewEuler(m1d.CFL, m1d.FinalTime, m1d.XMax, m1d.N, m1d.K, Euler1D.DFR_Roe, m1d.Case)
	case M_1DEulerDFR_LF:
		C = Euler1D.NewEuler(m1d.CFL, m1d.FinalTime, m1d.XMax, m1d.N, m1d.K, Euler1D.DFR_LaxFriedrichs, m1d.Case)
	case M_1DEulerDFR_Ave:
		C = Euler1D.NewEuler(m1d.CFL, m1d.FinalTime, m1d.XMax, m1d.N, m1d.K, Euler1D.DFR_Average, m1d.Case)
	case M_1DEuler:
		fallthrough
	default:
		C = Euler1D.NewEuler(m1d.CFL, m1d.FinalTime, m1d.XMax, m1d.N, m1d.K, Euler1D.Galerkin_LF, m1d.Case)
	}
	C.Run(m1d.Graph, m1d.Delay*time.Millisecond)
}

func LimitCFL(model ModelType1D, CFL float64) (CFLNew float64) {
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

func Defaults(model ModelType1D) (CFL, XMax float64, N, K, Case int) {
	return def_CFL[model], def_XMAX[model], def_N[model], def_K[model], def_CASE[model]
}

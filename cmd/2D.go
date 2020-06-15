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

import "C"
import (
	"fmt"
	"time"

	"github.com/notargets/gocfd/DG2D"
	"github.com/notargets/gocfd/model_problems/Euler1D"

	"github.com/spf13/cobra"
)

// TwoDCmd represents the 2D command
var TwoDCmd = &cobra.Command{
	Use:   "2D",
	Short: "Two dimensional solver, able to read grid files and output solutions",
	Long:  `Two dimensional solver, able to read grid files and output solutions`,
	Run: func(cmd *cobra.Command, args []string) {
		var (
			err error
		)
		fmt.Println("2D called")
		m2d := &Model2D{}
		mr, _ := cmd.Flags().GetInt("model")
		m2d.ModelRun = ModelType2D(mr)
		dr, _ := cmd.Flags().GetInt("delay")
		m2d.Delay = time.Duration(dr)
		Casep, _ := cmd.Flags().GetInt("case")
		m2d.Case = Euler1D.CaseType(Casep)
		m2d.FinalTime, _ = cmd.Flags().GetFloat64("finalTime")
		m2d.CFL, _ = cmd.Flags().GetFloat64("CFL")
		if m2d.GridFile, err = cmd.Flags().GetString("gridFile"); err != nil {
			panic(err)
		}
		m2d.Graph, _ = cmd.Flags().GetBool("graph")
		m2d.N, _ = cmd.Flags().GetInt("n")
		m2d.K, _ = cmd.Flags().GetInt("k")
		m2d.CFL = LimitCFL(ModelRun, CFL)
		Run2D(m2d)
	},
}

func init() {
	rootCmd.AddCommand(TwoDCmd)
	var CaseInt int
	CFL, XMax, N, K, CaseInt = Defaults(ModelRun)
	TwoDCmd.Flags().IntP("model", "m", int(ModelRun), "model to run: 0 = Advect1D, 1 = Maxwell1D, 2 = Euler1D")
	TwoDCmd.Flags().IntP("k", "k", K, "Number of elements in model")
	TwoDCmd.Flags().IntP("n", "n", N, "polynomial degree")
	TwoDCmd.Flags().IntP("delay", "d", 0, "milliseconds of delay for plotting")
	TwoDCmd.Flags().IntP("case", "c", int(CaseInt), "Case to run, for Euler: 0 = SOD Shock Tube, 1 = Density Wave")
	TwoDCmd.Flags().BoolP("graph", "g", false, "display a graph while computing solution")
	TwoDCmd.Flags().Float64("CFL", CFL, "CFL - increase for speedup, decrease for stability")
	TwoDCmd.Flags().Float64("finalTime", FinalTime, "FinalTime - the target end time for the sim")
	TwoDCmd.Flags().String("gridFile", "", "Grid file to read in Gambit (.neu) format")
}

type Model2D struct {
	K, N           int // Number of elements, Polynomial Degree
	Delay          time.Duration
	ModelRun       ModelType2D
	CFL, FinalTime float64
	Case           Euler1D.CaseType
	GridFile       string
	Graph          bool
}

type ModelType2D uint8

const (
	M_2DAdvect ModelType2D = iota
	M_2DMaxwell
	M_2DEuler
	M_2DAdvectDFR
	M_2DMaxwellDFR
	M_2DEulerDFR_Roe
	M_2DEulerDFR_LF
	M_2DEulerDFR_Ave
)

func Run2D(m2d *Model2D) {
	if len(m2d.GridFile) != 0 {
		DG2D.ReadGambit2d(m2d.GridFile, m2d.Graph)
	}
}

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
	"time"

	"github.com/notargets/gocfd/DG2D"

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
		m2d.FinalTime, _ = cmd.Flags().GetFloat64("finalTime")
		m2d.CFL, _ = cmd.Flags().GetFloat64("CFL")
		if m2d.GridFile, err = cmd.Flags().GetString("gridFile"); err != nil {
			panic(err)
		}
		m2d.Graph, _ = cmd.Flags().GetBool("graph")
		m2d.N, _ = cmd.Flags().GetInt("n")
		Run2D(m2d)
	},
}

func init() {
	rootCmd.AddCommand(TwoDCmd)
	var (
		CFL, FinalTime float64
		N              int
	)
	TwoDCmd.Flags().IntP("model", "m", int(0), "model to run")
	TwoDCmd.Flags().IntP("n", "n", N, "polynomial degree")
	TwoDCmd.Flags().IntP("delay", "d", 0, "milliseconds of delay for plotting")
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
	GridFile       string
	Graph          bool
}

type ModelType2D uint8

const (
	M_2DRoe ModelType2D = iota
)

func Run2D(m2d *Model2D) {
	_ = DG2D.NewNDG2D(m2d.N, m2d.GridFile, m2d.Graph)
}

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
	"io/ioutil"

	"github.com/notargets/gocfd/InputParameters"

	"github.com/notargets/gocfd/model_problems/Euler2D"

	"github.com/spf13/cobra"
)

type Model2D struct {
	GridFile          string
	NOrder            int
	ICFile            string
	ParallelProcLimit int
	Profile           bool
	fminP, fmaxP      *float64
}

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
		if m2d.GridFile, err = cmd.Flags().GetString("gridFile"); err != nil {
			panic(err)
		}
		m2d.NOrder, _ = cmd.Flags().GetInt("nOrder")
		if m2d.ICFile, err = cmd.Flags().GetString("inputConditionsFile"); err != nil {
			panic(err)
		}
		m2d.ParallelProcLimit, _ = cmd.Flags().GetInt("parallelProcs")
		m2d.Profile, _ = cmd.Flags().GetBool("profile")
		ip := processInput(m2d)
		Run2D(m2d, ip)
	},
}

func processInput(m2d *Model2D) (ip *InputParameters.InputParameters2D) {
	var (
		err error
	)
	if len(m2d.GridFile) == 0 {
		err := fmt.Errorf("must supply a grid file (-F, --gridFile) in .neu (Gambit neutral file) format")
		fmt.Printf("error: %s\n", err.Error())
	}
	if len(m2d.ICFile) == 0 {
		err := fmt.Errorf("must supply an input parameters file (-I, --inputConditionsFile) in .neu (Gambit neutral file) format")
		fmt.Printf("error: %s\n", err.Error())
		exampleFile := `
########################################
Title: "Test Case"
CFL: 1.
FluxType: Lax
InitType: IVortex # Can be "Freestream"
PolynomialOrder: 1
FinalTime: 4
########################################
`
		fmt.Printf("Example File Contents:%s\n", exampleFile)
	}
	ip = &InputParameters.InputParameters2D{}
	ip.Gamma = 1.4 // Default
	ip.Minf = 0.1  // Default
	var data []byte
	if len(m2d.ICFile) != 0 {
		if data, err = ioutil.ReadFile(m2d.ICFile); err != nil {
			panic(err)
		}
		if err = ip.Parse(data); err != nil {
			panic(err)
		}
	} else {
		ip.PolynomialOrder = m2d.NOrder // Default
		ip.FluxType = "Lax"             // Default
		ip.InitType = "Freestream"
	}
	return
}

func init() {
	rootCmd.AddCommand(TwoDCmd)
	TwoDCmd.Flags().StringP("gridFile", "F", "", "Grid file to read in Gambit (.neu) or SU2 (.su2) format")
	TwoDCmd.Flags().IntP("nOrder", "N", -1, "order of the polynomial used for the mesh")
	TwoDCmd.Flags().StringP("inputConditionsFile", "I", "", "YAML file for input parameters like:\n\t- CFL\n\t- NPR (nozzle pressure ratio)")
	TwoDCmd.Flags().IntP("parallelProcs", "p", 0, "limits the parallelism to the number of specified processes")
	TwoDCmd.Flags().Bool("profile", false, "generate a runtime profile of the solver, can be converted to PDF using 'go tool pprof -pdf filename'")
}

func Run2D(m2d *Model2D, ip *InputParameters.InputParameters2D) {
	c := Euler2D.NewEuler(ip, m2d.GridFile, m2d.ParallelProcLimit, true, m2d.Profile)
	c.SaveOutputMesh("meshfile.gcfd")
	c.Solve()
}

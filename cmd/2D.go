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
	"os"
	"sort"
	"time"

	"github.com/ghodss/yaml"

	"github.com/notargets/avs/chart2d"

	"github.com/notargets/gocfd/model_problems/Euler2D"

	"github.com/spf13/cobra"
)

type Model2D struct {
	GridFile   string
	ICFile     string
	Graph      bool
	GraphField int
	PlotSteps  int
	Delay      time.Duration
}

type InputParameters struct {
	Title           string                                `yaml:"Title"`
	CFL             float64                               `yaml:"CFL"`
	FluxType        string                                `yaml:"FluxType"`
	InitType        string                                `yaml:"InitType"`
	PolynomialOrder int                                   `yaml:"PolynomialOrder"`
	FinalTime       float64                               `yaml:"FinalTime"`
	BCs             map[string]map[int]map[string]float64 `yaml:"BCs"` // First key is BC name/type, second is parameter name
}

func (ip *InputParameters) Parse(data []byte) error {
	return yaml.Unmarshal(data, ip)
}

func (ip *InputParameters) Print() {
	fmt.Printf("\"%s\"\t\t= Title\n", ip.Title)
	fmt.Printf("%8.5f\t\t= CFL\n", ip.CFL)
	fmt.Printf("%8.5f\t\t= FinalTime\n", ip.FinalTime)
	fmt.Printf("[%s]\t\t\t= Flux Type\n", ip.FluxType)
	fmt.Printf("[%s]\t= InitType\n", ip.InitType)
	fmt.Printf("[%d]\t\t\t\t= Polynomial Order\n", ip.PolynomialOrder)
	keys := make([]string, len(ip.BCs))
	i := 0
	for k := range ip.BCs {
		keys[i] = k
		i++
	}
	sort.Strings(keys)
	for _, key := range keys {
		fmt.Printf("BCs[%s] = %v\n", key, ip.BCs[key])
	}
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
		if m2d.ICFile, err = cmd.Flags().GetString("inputConditionsFile"); err != nil {
			panic(err)
		}
		m2d.Graph, _ = cmd.Flags().GetBool("graph")
		m2d.GraphField, _ = cmd.Flags().GetInt("graphField")
		ps, _ := cmd.Flags().GetInt("plotSteps")
		m2d.PlotSteps = ps
		dr, _ := cmd.Flags().GetInt("delay")
		m2d.Delay = time.Duration(time.Duration(dr) * time.Millisecond)
		ip := processInput(m2d)
		Run2D(m2d, ip)
	},
}

func processInput(m2d *Model2D) (ip *InputParameters) {
	var (
		err      error
		willExit bool
	)
	if len(m2d.GridFile) == 0 {
		err := fmt.Errorf("must supply a grid file (-F, --gridFile) in .neu (Gambit neutral file) format")
		fmt.Printf("error: %s\n", err.Error())
		willExit = true
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
		fmt.Printf("Example File:%s\n", exampleFile)
		willExit = true
	}
	if willExit {
		os.Exit(1)
	}
	if len(m2d.ICFile) != 0 {
		var data []byte
		if data, err = ioutil.ReadFile(m2d.ICFile); err != nil {
			panic(err)
		}
		ip = &InputParameters{}
		if err = ip.Parse(data); err != nil {
			panic(err)
		}
	}
	return
}

func init() {
	rootCmd.AddCommand(TwoDCmd)
	TwoDCmd.Flags().StringP("gridFile", "F", "", "Grid file to read in Gambit (.neu) format")
	TwoDCmd.Flags().StringP("inputConditionsFile", "I", "", "YAML file for input parameters like:\n\t- CFL\n\t- NPR (nozzle pressure ratio)")
	TwoDCmd.Flags().BoolP("graph", "g", false, "display a graph while computing solution")
	TwoDCmd.Flags().IntP("delay", "d", 0, "milliseconds of delay for plotting")
	TwoDCmd.Flags().IntP("plotSteps", "s", 1, "number of steps before plotting each frame")
	TwoDCmd.Flags().IntP("graphField", "q", 0, "which field should be displayed - 0=density, 1,2=momenta, 3=energy")
}

func Run2D(m2d *Model2D, ip *InputParameters) {
	c := Euler2D.NewEuler(
		ip.FinalTime, ip.PolynomialOrder, m2d.GridFile, ip.CFL,
		Euler2D.NewFluxType(ip.FluxType), Euler2D.NewInitType(ip.InitType),
		false, true)
	//m2d.FinalTime, m2d.N, m2d.GridFile, m2d.CFL, m2d.FluxType, m2d.InitType,
	pm := &Euler2D.PlotMeta{
		Plot:            m2d.Graph,
		Scale:           1.1,
		Field:           Euler2D.PlotField(m2d.GraphField),
		FieldMinP:       nil,
		FieldMaxP:       nil,
		FrameTime:       m2d.Delay,
		StepsBeforePlot: m2d.PlotSteps,
		LineType:        chart2d.NoLine,
	}
	c.Solve(pm)
}

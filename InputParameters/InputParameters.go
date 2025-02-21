package InputParameters

import (
	"fmt"
	"sort"

	"github.com/ghodss/yaml"
)

// Parameters obtained from the YAML input file
type InputParameters2D struct {
	Title             string                                `yaml:"Title"`
	CFL               float64                               `yaml:"CFL"`
	FluxType          string                                `yaml:"FluxType"`
	InitType          string                                `yaml:"InitType"`
	PolynomialOrder   int                                   `yaml:"PolynomialOrder"`
	FinalTime         float64                               `yaml:"FinalTime"`
	Minf              float64                               `yaml:"Minf"`
	Gamma             float64                               `yaml:"Gamma"`
	Alpha             float64                               `yaml:"Alpha"`
	BCs               map[string]map[int]map[string]float64 `yaml:"BCs"` // First key is BC name/type, second is parameter name
	LocalTimeStepping bool                                  `yaml:"LocalTimeStep"`
	MaxIterations     int                                   `yaml:"MaxIterations"`
	ImplicitSolver    bool                                  `yaml:"ImplicitSolver"`
	Limiter           string                                `yaml:"Limiter"`
	Kappa             float64                               `yaml:"Kappa"`
}

func (ip *InputParameters2D) Parse(data []byte) error {
	return yaml.Unmarshal(data, ip)
}

func (ip *InputParameters2D) Print() {
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

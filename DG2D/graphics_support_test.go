package DG2D

import (
	"fmt"
	"testing"

	"github.com/notargets/gocfd/utils"

	"github.com/stretchr/testify/assert"
)

func TestPlotVariousFields(t *testing.T) {
	var (
		N = 3
	)
	if !testing.Verbose() {
		return
	}
	// angle := 82.8
	// dfr := CreateEquiTriMesh(N, angle)
	// gm := dfr.CreateAVSGraphMesh()
	// dfr := NewDFR2D(N, false, "test_data/test_tris_9.neu")
	dfr := NewDFR2D(N, false, "test_data/test_10tris_centered.neu")

	// SolutionX and S are dimension (Np,K)
	// Np := dfr.SolutionElement.Np
	// K := dfr.K
	// for k := 0; k < K; k++ {
	// 	dfr.SolutionX.Col(k).Transpose().String("XElement:" + strconv.Itoa(k))
	// 	dfr.SolutionY.Col(k).Transpose().String("YElement:" + strconv.Itoa(k))
	// 	minX := dfr.SolutionX.Col(k).Min()
	// 	maxX := dfr.SolutionX.Col(k).Max()
	// 	if minX < 0 && maxX < 0 {
	// 		fmt.Printf("Element [%d] is pre-shock\n\n", k)
	// 	} else if minX < 0 && maxX >= 0 {
	// 		fmt.Printf("Element [%d] is mixed pre- and POST-shock\n\n", k)
	// 	} else {
	// 		fmt.Printf("Element [%d] is POST-shock\n\n", k)
	// 	}
	// }
	X, Y := dfr.SolutionX.DataP, dfr.SolutionY.DataP
	// field := SetTestField(R, S, RADIAL2TEST)
	// field := SetTestField(R, S, FIXEDVORTEXTEST)
	// field := SetTestField(R, S, NORMALSHOCKTESTM12)
	// field := SetTestField(R, S, NORMALSHOCKTESTM2)
	field := SetTestField(X, Y, NORMALSHOCKTESTM5)
	fM := utils.NewMatrix(dfr.SolutionElement.Np,
		dfr.K, field[:dfr.K*dfr.SolutionElement.Np])
	fM.PrintDims("fM Before Mod")
	// fM.String("fM Before Mod Filter")
	// dfr.FilterMod.PrintDims("FilterMod")
	// fM = fM.Transpose().Mul(dfr.FilterMod).Transpose()
	// fM.PrintDims("After Mod Filter")
	// fM.String("After Mod Filter")
	fMin, fMax := getFieldMinMax(field[:dfr.K*dfr.SolutionElement.Np])
	// dfr.GraphInterp.PrintDims("GraphInterp")
	// The graph mesh field dimensions are (K,Np)
	fMGraph := dfr.GraphInterp.Mul(fM).Transpose()
	// fMGraph := dfr.GraphInterpMod.Mul(fM).Transpose()
	fMGraph.PrintDims("Graph Field Dims")
	fMGraph.Print("Graph Field")
	fmt.Printf("fMin: %f, fMax: %f\n", fMin, fMax)
	PlotField(fMGraph.DataP, dfr.GraphMesh, fMin, fMax)
}

func TestDFR2D_WriteAVSGraphMesh(t *testing.T) {
	if !testing.Verbose() {
		return
	}
	dfr := NewDFR2D(3, false, "test_data/test_tris_9.neu")
	dfr.OutputMesh("test_data/test_tris_9.gobcfd", nil)
	md, gm, BCXY, err := ReadMesh("test_data/test_tris_9.gobcfd")
	_, _, _ = md, gm, BCXY
	assert.NoError(t, err)
	assert.NotNil(t, md)
	assert.Equal(t, 10, md.NumBaseElements)
	assert.Equal(t, 36, md.NumPerElement)
	assert.Equal(t, 360, md.NumElements)
	assert.Equal(t, 3, md.Order)
	fmt.Println(md)
	// fmt.Printf("Git tag: %v\n", md.GitVersion)
	// md, gm, BCXY, err := ReadMesh("test_data/meshfile.gobcfd")
	// for name, XYranges := range BCXY {
	// 	fmt.Println(name)
	// 	for i, xy := range XYranges {
	// 		fmt.Println("I range: ", i, "range: ", xy)
	// 	}
	// }
	// PlotTriMesh(gm)
	// for {
	// }
}

func TestCreateAVSGraphMesh(t *testing.T) {
	var (
		N = 7
	)
	if !testing.Verbose() {
		return
	}
	// dfr := NewDFR2D(N, false, "test_data/test_tris_5.neu")
	// dfr := NewDFR2D(N, false, "test_data/test_tris_1tri.neu")
	dfr := NewDFR2D(N, false, "test_data/test_tris_9.neu")
	gm := dfr.CreateAVSGraphMesh()
	PlotTriMesh(gm)
}

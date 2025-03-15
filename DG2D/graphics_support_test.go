package DG2D

import (
	"fmt"
	"math"
	"testing"

	"github.com/notargets/gocfd/utils"

	"github.com/stretchr/testify/assert"

	"github.com/notargets/avs/chart2d"
	"github.com/notargets/avs/geometry"
	utils2 "github.com/notargets/avs/utils"
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

	// SolutionX and Y are dimension (Np,K)
	// Np := dfr.SolutionElement.Np
	// K := dfr.K
	// for k := 0; k < K; k++ {
	// 	dfr.SolutionX.Col(k).Transpose().Print("XElement:" + strconv.Itoa(k))
	// 	dfr.SolutionY.Col(k).Transpose().Print("YElement:" + strconv.Itoa(k))
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
	// field := setTestField(X, Y, RADIAL2TEST)
	// field := setTestField(X, Y, FIXEDVORTEXTEST)
	// field := setTestField(X, Y, NORMALSHOCKTESTM12)
	// field := setTestField(X, Y, NORMALSHOCKTESTM2)
	field := setTestField(X, Y, NORMALSHOCKTESTM5)
	fM := utils.NewMatrix(dfr.SolutionElement.Np,
		dfr.K, field[:dfr.K*dfr.SolutionElement.Np])
	fM.PrintDims("fM Before Mod")
	// fM.Print("fM Before Mod Filter")
	// dfr.FilterMod.PrintDims("FilterMod")
	// fM = fM.Transpose().Mul(dfr.FilterMod).Transpose()
	// fM.PrintDims("After Mod Filter")
	// fM.Print("After Mod Filter")
	fMin, fMax := getFieldMinMax(field[:dfr.K*dfr.SolutionElement.Np])
	// dfr.GraphInterp.PrintDims("GraphInterp")
	// The graph mesh field dimensions are (K,Np)
	fMGraph := dfr.GraphInterp.Mul(fM).Transpose()
	// fMGraph := dfr.GraphInterpMod.Mul(fM).Transpose()
	fMGraph.PrintDims("Graph Field Dims")
	fMGraph.Print("Graph Field")
	fmt.Printf("fMin: %f, fMax: %f\n", fMin, fMax)
	plotField(fMGraph.DataP, dfr.GraphMesh, fMin, fMax)
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
	// plotMesh(gm)
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
	plotMesh(gm)
}

func plotMesh(gm geometry.TriMesh) {
	xMin, xMax, yMin, yMax := getMinMax(gm)
	ch := chart2d.NewChart2D(xMin, xMax, yMin, yMax,
		1024, 1024, utils2.WHITE, utils2.BLACK)
	// Create a vector field including the three vertices
	ch.AddTriMesh(gm)
	for {
	}
}

func getMinMax(gm geometry.TriMesh) (xMin, xMax, yMin, yMax float32) {
	var (
		x, y  float32
		lenXY = len(gm.XY) / 2
	)
	for i := 0; i < lenXY; i++ {
		x, y = gm.XY[i*2+0], gm.XY[i*2+1]
		if i == 0 {
			xMin = x
			xMax = x
			yMin = y
			yMax = y
		} else {
			if x < xMin {
				xMin = x
			}
			if x > xMax {
				xMax = x
			}
			if y < yMin {
				yMin = y
			}
			if y > yMax {
				yMax = y
			}
		}
	}
	return
}

func getFieldMinMax(field []float64) (fMin, fMax float64) {
	for i, f := range field {
		if i == 0 {
			fMin = f
			fMax = f
		}
		if f < fMin {
			fMin = f
		}
		if f > fMax {
			fMax = f
		}
	}
	return
}

func plotField(field []float64, gm geometry.TriMesh, FMin, FMax float64,
	xMM ...float64) {
	var xMin, xMax, yMin, yMax float32

	if len(xMM) == 4 {
		xMin, xMax = float32(xMM[0]), float32(xMM[1])
		yMin, yMax = float32(xMM[2]), float32(xMM[3])
	} else {
		xMin, xMax = -1, 1
		yMin, yMax = -1, 1
	}
	ch := chart2d.NewChart2D(xMin, xMax, yMin, yMax,
		1024, 1024, utils2.WHITE, utils2.BLACK)
	// Create a vector field including the three vertices
	var pField []float32
	var fMin, fMax float32
	fMin, fMax = math.MaxFloat32, -math.MaxFloat32
	pField = make([]float32, len(field))
	for i, f := range field {
		f32 := float32(f)
		if fMin > f32 {
			fMin = f32
		}
		if fMax < f32 {
			fMax = f32
		}
		pField[i] = float32(f)
	}
	vs := geometry.VertexScalar{
		TMesh:       &gm,
		FieldValues: pField,
	}
	fmt.Printf("Interpolated fMin: %f, fMax: %f\n", fMin, fMax)
	ch.AddShadedVertexScalar(&vs, float32(FMin), float32(FMax))
	ch.AddTriMesh(gm)
	line := []float32{0, -5, 0, 5, -5, 0, 5, 0}
	ch.AddLine(line, utils2.RED)
	for {
	}
}

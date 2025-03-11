package DG2D

import (
	"fmt"
	"math"
	"testing"

	"github.com/stretchr/testify/assert"

	"github.com/notargets/avs/chart2d"
	"github.com/notargets/avs/geometry"
	utils2 "github.com/notargets/avs/utils"
)

func plotField(field []float64, fieldNum int, gm geometry.TriMesh,
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
	vertsLen := len(field) / 4
	var fMin, fMax float32
	fMin, fMax = math.MaxFloat32, -math.MaxFloat32
	n := 0 // Density
	n = 2  // V momentum
	n = 3  // rho*E energy
	n = 0  // Density
	n = 1  // U momentum
	n = 3  // rho*E energy
	n = 0  // Density
	for i := 0; i < vertsLen; i++ {
		pField = append(pField, float32(field[i+n*vertsLen]))
		// fmt.Printf("F[%d] = %.2f\n", i, pField[i])
		if pField[i] < fMin {
			fMin = pField[i]
		}
		if pField[i] > fMax {
			fMax = pField[i]
		}
	}
	vs := geometry.VertexScalar{
		TMesh:       &gm,
		FieldValues: pField,
	}
	_ = vs
	_ = ch
	// ch.AddContourVertexScalar(&vs, fMin, fMax, 100)
	ch.AddShadedVertexScalar(&vs, fMin, fMax)
	ch.AddTriMesh(gm)
	for {
	}
}

func TestPlotVariousFields(t *testing.T) {
	var (
		N = 3
	)
	if !testing.Verbose() {
		return
	}
	angle := 82.8
	dfr := CreateEquiTriMesh(N, angle)
	gm := dfr.CreateAVSGraphMesh()
	X, Y := convXYtoXandY(gm.XY)
	field := setTestField(X, Y, NORMALSHOCKTESTM5)
	// field := setTestField(X, Y, NORMALSHOCKTESTM2)
	// field := setTestField(X, Y, RADIAL2TEST)
	n := 0 // Density
	n = 1  // U momentum
	n = 2  // V momentum
	n = 3  // rho*E energy
	n = 0  // Density
	plotField(field, n, gm)
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

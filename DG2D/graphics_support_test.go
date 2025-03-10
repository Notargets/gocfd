package DG2D

import (
	"fmt"
	"testing"

	"github.com/stretchr/testify/assert"

	"github.com/notargets/avs/chart2d"
	"github.com/notargets/avs/geometry"
	utils2 "github.com/notargets/avs/utils"
)

func TestDFR2D_WriteAVSGraphMesh(t *testing.T) {
	if !testing.Verbose() {
		return
	}
	dfr := NewDFR2D(3, false, "test_data/test_tris_9.neu")
	dfr.WriteAVSGraphMesh(dfr.CreateAVSGraphMesh(), "test_data/test_tris_9.gobcfd")
	md, gm, err := ReadMesh("test_data/test_tris_9.gobcfd")
	_, _ = md, gm
	assert.NoError(t, err)
	fmt.Println(md)
	fmt.Printf("Git tag: %v\n", md.GitVersion)
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

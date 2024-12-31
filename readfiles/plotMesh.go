package readfiles

import (
	"fmt"
	"image/color"

	"github.com/notargets/gocfd/InputParameters"

	"github.com/notargets/avs/chart2d"
	graphics2D "github.com/notargets/avs/geometry"
	utils2 "github.com/notargets/avs/utils"
	"github.com/notargets/gocfd/types"
	"github.com/notargets/gocfd/utils"
)

func PrintTriMesh(trimesh *graphics2D.TriMesh) {
	fmt.Println("Geometry, X, Y for each vertex:")
	for _, pt := range trimesh.Geometry {
		fmt.Printf("%5.2f %5.2f\n", pt.X[0], pt.X[1])
	}
	fmt.Println("Vertex indices, 3 for each triangle:")
	for _, t := range trimesh.Triangles {
		for i := 0; i < 3; i++ {
			fmt.Printf("%d ", t.Nodes[i])
		}
		fmt.Printf("\n")
	}
	fmt.Println("Scalar color value, one for each index in the triangle:")
	for _, a := range trimesh.Attributes {
		for i := 0; i < 3; i++ {
			fmt.Printf("%5.2f ", a[i])
		}
		fmt.Printf("\n")
	}
}

func PlotMesh(VX, VY utils.Vector, EToV, X, Y utils.Matrix, plotPoints bool, pm *InputParameters.PlotMeta) (chart *chart2d.Chart2D) {
	var (
		points   []graphics2D.Point
		trimesh  graphics2D.TriMesh
		vxD, vyD = VX.DataP, VY.DataP
		K, _     = EToV.Dims()
	)
	points = make([]graphics2D.Point, VX.Len())
	for i, vx := range vxD {
		points[i].X[0] = float32(vx)
		points[i].X[1] = float32(vyD[i])
	}
	trimesh.Triangles = make([]graphics2D.Triangle, K)
	colorMap := utils2.NewColorMap(0, float32(types.BC_PeriodicReversed), 1)
	trimesh.Attributes = make([][]float32, K) // One BC attribute per face
	for k := 0; k < K; k++ {
		trimesh.Attributes[k] = make([]float32, 3)
		for i := 0; i < 3; i++ {
			trimesh.Triangles[k].Nodes[i] = int32(EToV.At(k, i))
			trimesh.Attributes[k][i] = float32(types.BC_None)
		}
	}
	trimesh.Geometry = points

	PrintTriMesh(&trimesh)

	box := graphics2D.NewBoundingBox(trimesh.GetGeometry())
	box = box.Translate([2]float32{float32(pm.TranslateX), float32(pm.TranslateY)}).Scale(float32(pm.Scale))
	chart = chart2d.NewChart2D(1920, 1920, box.XMin[0], box.XMax[0], box.XMin[1], box.XMax[1])
	chart.AddColorMap(colorMap)
	go chart.Plot()
	white := color.RGBA{
		R: 255,
		G: 255,
		B: 255,
		A: 0,
	}
	//black := color.RGBA{
	//	R: 0,
	//	G: 0,
	//	B: 0,
	//	A: 0,
	//}
	_, _ = colorMap, white
	if err := chart.AddTriMesh("TriMesh", trimesh,
		chart2d.NoGlyph, 0.005, chart2d.Solid, white); err != nil {
		panic("unable to add graph series")
	}
	//var ptsGlyph chart2d.GlyphType
	//ptsGlyph = chart2d.NoGlyph
	//if plotPoints {
	//	ptsGlyph = chart2d.CircleGlyph
	//}
	//if err := chart.AddSeries("Elements", X.Transpose().DataP, Y.Transpose().DataP,
	//	ptsGlyph, 0.0005, chart2d.NoLine, black); err != nil {
	//	panic(err)
	//}

	return
}

func PlotMesh_old(VX, VY utils.Vector, EToV, X, Y utils.Matrix, plotPoints bool, pm *InputParameters.PlotMeta) (chart *chart2d.Chart2D) {
	var (
		points   []graphics2D.Point
		trimesh  graphics2D.TriMesh
		vxD, vyD = VX.DataP, VY.DataP
		K, _     = EToV.Dims()
	)
	points = make([]graphics2D.Point, VX.Len())
	for i, vx := range vxD {
		points[i].X[0] = float32(vx)
		points[i].X[1] = float32(vyD[i])
	}
	trimesh.Triangles = make([]graphics2D.Triangle, K)
	colorMap := utils2.NewColorMap(0, float32(types.BC_PeriodicReversed), 1)
	trimesh.Attributes = make([][]float32, K) // One BC attribute per face
	for k := 0; k < K; k++ {
		trimesh.Attributes[k] = make([]float32, 3)
		for i := 0; i < 3; i++ {
			trimesh.Triangles[k].Nodes[i] = int32(EToV.At(k, i))
			trimesh.Attributes[k][i] = float32(types.BC_None)
		}
	}
	trimesh.Geometry = points
	box := graphics2D.NewBoundingBox(trimesh.GetGeometry())
	box = box.Translate([2]float32{float32(pm.TranslateX), float32(pm.TranslateY)}).Scale(float32(pm.Scale))
	chart = chart2d.NewChart2D(1920, 1920, box.XMin[0], box.XMax[0], box.XMin[1], box.XMax[1])
	chart.AddColorMap(colorMap)
	go chart.Plot()
	white := color.RGBA{
		R: 255,
		G: 255,
		B: 255,
		A: 0,
	}
	black := color.RGBA{
		R: 0,
		G: 0,
		B: 0,
		A: 0,
	}
	_, _ = colorMap, white
	if err := chart.AddTriMesh("TriMesh", trimesh,
		chart2d.NoGlyph, 0.005, chart2d.Solid, white); err != nil {
		panic("unable to add graph series")
	}
	var ptsGlyph chart2d.GlyphType
	ptsGlyph = chart2d.NoGlyph
	if plotPoints {
		ptsGlyph = chart2d.CircleGlyph
	}
	if err := chart.AddSeries("Elements", X.Transpose().DataP, Y.Transpose().DataP,
		ptsGlyph, 0.0005, chart2d.NoLine, black); err != nil {
		panic(err)
	}

	return
}

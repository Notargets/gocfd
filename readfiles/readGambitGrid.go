package readfiles

import (
	"bufio"
	"fmt"
	"image/color"
	"io"
	"os"
	"strings"

	"github.com/notargets/gocfd/types"

	"github.com/notargets/avs/chart2d"
	utils2 "github.com/notargets/avs/utils"

	graphics2D "github.com/notargets/avs/geometry"

	"github.com/notargets/gocfd/utils"
)

type Material struct {
	ElementCount  int
	MaterialValue float64
	Title         string
}

func ReadGambit2d(filename string, verbose bool) (K int, VX, VY utils.Vector, EToV utils.Matrix, BCEdges types.BCMAP) {
	var (
		file   *os.File
		err    error
		reader *bufio.Reader
	)
	if verbose {
		fmt.Printf("Reading Gambit Neutral file named: %s\n", filename)
	}
	if file, err = os.Open(filename); err != nil {
		panic(fmt.Errorf("unable to open file %s\n %s", filename, err))
	}
	defer file.Close()
	reader = bufio.NewReader(file)

	// Skip firest six lines
	skipLines(6, reader)

	// Get dimensions
	Nv, K, Nmats, Nbcs, Nsd := ReadHeader(reader)
	skipLines(2, reader)

	if verbose {
		fmt.Printf("Nv = %d, K = %d\n", Nv, K)
		fmt.Printf("Nmats = %d, Nbcs = %d\n%d space dimensions\n", Nmats, Nbcs, Nsd)
	}
	if Nsd > 3 || Nsd < 2 {
		panic("space dimensions not 2 or 3")
	}

	NFaces, bIs3D, bCoord3D, bElement3D, bTET := CalculateGeomAttributes(Nsd, false)

	if verbose {
		fmt.Printf("NFaces = %d, bIs3D is %v, bCoord3D is %v, bElement3D is %v, bTET is %v\n",
			NFaces, bIs3D, bCoord3D, bElement3D, bTET)
	}

	var VZ utils.Vector
	if bCoord3D {
		VX, VY, VZ = Read3DVertices(Nv, reader)
	} else {
		VX, VY = Read2DVertices(Nv, reader)
	}
	skipLines(2, reader)

	// Read Elements
	if bTET {
		EToV = ReadTets(K, reader)
	} else {
		EToV = ReadTris(K, reader)
	}
	skipLines(2, reader)

	if verbose {
		switch Nsd {
		case 2:
			fmt.Printf("Bounding Box:\nXMin/XMax = %5.3f, %5.3f\nYMin/YMax = %5.3f, %5.3f\n",
				VX.Min(), VX.Max(), VY.Min(), VY.Max())
		case 3:
			fmt.Printf("Bounding Box:\nXMin/XMax = %5.3f, %5.3f\nYMin/YMax = %5.3f, %5.3f\nZMin/ZMax = %5.3f, %5.3f\n",
				VX.Min(), VX.Max(), VY.Min(), VY.Max(), VZ.Min(), VZ.Max())
		}
	}

	matGroups := make(map[int]*Material)
	// Read material values
	epsilon := utils.NewVector(K)
	for i := 0; i < Nmats; i++ {
		gn, elnum, matval, title := ReadMaterialHeader(reader)
		matGroups[gn] = &Material{
			ElementCount:  elnum,
			MaterialValue: matval,
			Title:         title,
		}
		ReadMaterialGroup(reader, elnum, matval, epsilon, false)
		skipLines(2, reader)
	}

	// Read BCs
	BCEdges = ReadBCS(Nbcs, K, NFaces, reader, EToV)
	return
}

func PlotMesh(VX, VY utils.Vector, EToV, X, Y utils.Matrix, plotPoints bool) (chart *chart2d.Chart2D) {
	var (
		points   []graphics2D.Point
		trimesh  graphics2D.TriMesh
		vxD, vyD = VX.Data(), VY.Data()
		K, _     = EToV.Dims()
	)
	points = make([]graphics2D.Point, VX.Len())
	for i, vx := range vxD {
		points[i].X[0] = float32(vx)
		points[i].X[1] = float32(vyD[i])
	}
	trimesh.Triangles = make([]graphics2D.Triangle, K)
	colorMap := utils2.NewColorMap(0, float32(types.BC_Out), 1)
	trimesh.Attributes = make([][]float32, K) // One BC attribute per face
	for k := 0; k < K; k++ {
		trimesh.Triangles[k].Nodes[0] = int32(EToV.At(k, 0))
		trimesh.Triangles[k].Nodes[1] = int32(EToV.At(k, 1))
		trimesh.Triangles[k].Nodes[2] = int32(EToV.At(k, 2))
	}
	trimesh.Geometry = points
	box := graphics2D.NewBoundingBox(trimesh.GetGeometry())
	box = box.Scale(1.5)
	chart = chart2d.NewChart2D(1920, 1920, box.XMin[0], box.XMax[0], box.XMin[1], box.XMax[1])
	chart.AddColorMap(colorMap)
	go chart.Plot()
	white := color.RGBA{
		R: 255,
		G: 255,
		B: 255,
		A: 0,
	}
	_ = colorMap
	if err := chart.AddTriMesh("TriMesh", trimesh,
		chart2d.CrossGlyph, chart2d.Solid, white); err != nil {
		panic("unable to add graph series")
	}
	var ptsGlyph chart2d.GlyphType
	ptsGlyph = chart2d.NoGlyph
	if plotPoints {
		ptsGlyph = chart2d.CircleGlyph
	}
	if err := chart.AddSeries("Elements", X.Transpose().Data(), Y.Transpose().Data(),
		ptsGlyph, chart2d.NoLine, white); err != nil {
		panic(err)
	}

	return
}

func ReadBCS(Nbcs, K, NFaces int, reader *bufio.Reader, EToV utils.Matrix) (BCEdges types.BCMAP) {
	var (
		line, bctyp string
		err         error
		nargs       int
		n, bcid     int
	)
	BCEdges = make(types.BCMAP, Nbcs)
	for i := 0; i < Nbcs; i++ {
		// Read BC header, if BC text is "Cyl", read a float parameter
		if i != 0 {
			skipLines(1, reader)
		}
		line = getLine(reader)
		if n, err = fmt.Sscanf(line, "%32s", &bctyp); err != nil {
			panic(err)
		}
		bctyp = strings.ToLower(strings.Trim(bctyp, " "))
		var paramf float64
		var numfaces int
		switch bctyp {
		case "cyl":
			if n, err = fmt.Sscanf(line, "%32s%8f%8d", &bctyp, &paramf, &numfaces); err != nil {
				panic(err)
			}
		default:
			if n, err = fmt.Sscanf(line, "%32s%8d%8d", &bctyp, &bcid, &numfaces); err != nil {
				panic(err)
			}
		}
		edges := make([]types.EdgeInt, numfaces)
		var verts [3]int
		for i := 0; i < numfaces; i++ {
			line = getLine(reader)
			nargs = 3
			var kp1, n2, faceNumberp1 int
			if n, err = fmt.Sscanf(line, "%d %d %d", &kp1, &n2, &faceNumberp1); err != nil || n < nargs {
				if err == nil && n < nargs {
					err = fmt.Errorf("read fewer than required dimensions, read %d, need %d\n, line: %s", n, nargs, line)
				}
				panic(err)
			}
			verts[0] = int(EToV.At(kp1-1, 0))
			verts[1] = int(EToV.At(kp1-1, 1))
			verts[2] = int(EToV.At(kp1-1, 2))
			var e types.EdgeInt
			switch faceNumberp1 {
			case 1:
				e = types.NewEdgeInt([2]int{verts[0], verts[1]})
			case 2:
				e = types.NewEdgeInt([2]int{verts[1], verts[2]})
			case 3:
				e = types.NewEdgeInt([2]int{verts[2], verts[0]})
			}
			edges[i] = e
		}
		BCEdges.AddEdges(types.NewBCTAG(bctyp), edges)
		skipLines(1, reader)
	}
	return
}

func ReadMaterialGroup(reader *bufio.Reader, elementCount int, matval float64, epsilon utils.Vector, verbose bool) {
	var (
		n       int
		nn      = make([]int, 10)
		epsData = epsilon.Data()
		err     error
		added   int
	)
	if elementCount%10 != 0 {
		added = 1
	}
	numLines := elementCount/10 + added
	if verbose {
		fmt.Printf("Reading %d lines of materials with %d elements\n", numLines, elementCount)
	}
	for i := 0; i < numLines; i++ {
		line := getLine(reader)
		nargs := 10
		if n, err = fmt.Sscanf(line, "%d %d %d %d %d %d %d %d %d %d", &nn[0], &nn[1], &nn[2], &nn[3], &nn[4], &nn[5], &nn[6], &nn[7], &nn[8], &nn[9]); err != nil || n < nargs {
			if !(n < nargs && i == numLines-1) {
				if err == nil && n < nargs {
					err = fmt.Errorf("read fewer than %d dimensions, read %d, line: %s", nargs, n, line)
				}
				panic(err)
			}
		}
		for j := 0; j < n; j++ {
			epsData[nn[j]-1] = matval
		}
	}
	return
}

func ReadMaterialHeader(reader *bufio.Reader) (gn, elnum int, matval float64, title string) {
	/*
	   GROUP:           1 ELEMENTS:        977 MATERIAL:      1.000 NFLAGS:          0
	                     epsilon: 1.000
	          0
	*/
	var (
		line = getLine(reader)
		n    int
		err  error
	)
	nargs := 3
	if n, err = fmt.Sscanf(line, "GROUP: %11d ELEMENTS:%11d MATERIAL:%11f", &gn, &elnum, &matval); err != nil || n < nargs {
		if err == nil && n < nargs {
			err = fmt.Errorf("read fewer than %d dimensions, read %d, line: %s", nargs, n, line)
		}
		panic(err)
	}
	title = getLine(reader)
	skipLines(1, reader)
	return
}

func ReadHeader(reader *bufio.Reader) (Nv, K, Nmats, Nbcs, Nsd int) {
	/*
		Nv      // num nodes in mesh
		K       // num elements
		Nmats   // num material groups
		Nbcs    // num boundary groups
		Nsd;    // num space dimensions
	*/
	var (
		line   = getLine(reader)
		n, dum int
		err    error
	)
	nargs := 6
	if n, err = fmt.Sscanf(line, "%d %d %d %d %d %d", &Nv, &K, &Nmats, &Nbcs, &Nsd, &dum); err != nil || n < nargs {
		if err == nil && n < nargs {
			err = fmt.Errorf("read fewer than %d dimensions, read %d, line: %s", nargs, n, line)
		}
		panic(err)
	}
	return
}

func Read2DVertices(Nv int, reader *bufio.Reader) (VX, VY utils.Vector) {
	var (
		line   string
		err    error
		n, ind int
	)
	nargs := 3
	VX, VY = utils.NewVector(Nv), utils.NewVector(Nv)
	vx, vy := VX.Data(), VY.Data()
	for i := 0; i < Nv; i++ {
		line = getLine(reader)
		if n, err = fmt.Sscanf(line, "%d", &ind); err != nil || n < 1 {
			if err == nil && n < nargs {
				err = fmt.Errorf("read fewer than required dimensions, read %d, need %d\n, line: %s", n, nargs, line)
			}
			err = fmt.Errorf("error reading index, line: %s, err: %s\n", line, err.Error())
			panic(err)
		}
		if n, err = fmt.Sscanf(line, "%d %f %f", &ind, &vx[ind-1], &vy[ind-1]); err != nil || n < nargs {
			if err == nil && n < nargs {
				err = fmt.Errorf("read fewer than required dimensions, read %d, need %d\n, line: %s", n, nargs, line)
			}
			panic(err)
		}
	}
	return
}

func Read3DVertices(Nv int, reader *bufio.Reader) (VX, VY, VZ utils.Vector) {
	var (
		line   string
		err    error
		n, ind int
	)
	nargs := 4
	VX, VY, VZ = utils.NewVector(Nv), utils.NewVector(Nv), utils.NewVector(Nv)
	vx, vy, vz := VX.Data(), VY.Data(), VZ.Data()
	for i := 0; i < Nv; i++ {
		line = getLine(reader)
		if n, err = fmt.Sscanf(line, "%d", &ind); err != nil || n < 1 {
			if err == nil && n < nargs {
				err = fmt.Errorf("read fewer than required dimensions, read %d, need %d\n, line: %s", n, nargs, line)
			}
			panic(err)
		}
		if n, err = fmt.Sscanf(line, "%d %f %f %f", &ind, &vx[ind-1], &vy[ind-1], &vz[ind-1]); err != nil || n < nargs {
			if err == nil && n < nargs {
				err = fmt.Errorf("read fewer than required dimensions, read %d, need %d\n, line: %s", n, nargs, line)
			}
			panic(err)
		}
	}
	return
}

func ReadTets(K int, reader *bufio.Reader) (EToV utils.Matrix) {
	//---------------------------------------------
	// Tetrahedra in 3D:
	//---------------------------------------------
	// ENDOFSECTION
	//    ELEMENTS/CELLS 1.3.0
	//     1  6  4      248     247     385     265
	//     2  6  4      248     249     273     397
	//---------------------------------------------
	var (
		line                       string
		err                        error
		n, ind, typ, nfaces, nargs int
	)
	EToV = utils.NewMatrix(K, 4)
	for i := 0; i < K; i++ {
		line = getLine(reader)
		nargs = 7
		var n1, n2, n3, n4 int
		if n, err = fmt.Sscanf(line, "%d %d %d %d %d %d %d", &ind, &typ, &nfaces, &n1, &n2, &n3, &n4); err != nil || n < nargs {
			if err == nil && n < nargs {
				err = fmt.Errorf("read fewer than required dimensions, read %d, need %d\n, line: %s", n, nargs, line)
			}
			panic(err)
		}
		EToV.Set(ind-1, 0, float64(n1))
		EToV.Set(ind-1, 1, float64(n2))
		EToV.Set(ind-1, 2, float64(n3))
		EToV.Set(ind-1, 3, float64(n4))
	}
	return
}

func ReadTris(K int, reader *bufio.Reader) (EToV utils.Matrix) {
	//-------------------------------------
	// Triangles in 3D:
	//-------------------------------------
	// ENDOFSECTION
	//    ELEMENTS/CELLS 1.3.0
	//      1  3  3        1       2       3
	//      2  3  3        3       2       4
	var (
		line                       string
		err                        error
		n, ind, typ, nfaces, nargs int
	)
	EToV = utils.NewMatrix(K, 3)
	for i := 0; i < K; i++ {
		line = getLine(reader)
		nargs = 6
		var n1, n2, n3 int
		if n, err = fmt.Sscanf(line, "%d %d %d %d %d %d", &ind, &typ, &nfaces, &n1, &n2, &n3); err != nil || n < nargs {
			if err == nil && n < nargs {
				err = fmt.Errorf("read fewer than required dimensions, read %d, need %d\n, line: %s", n, nargs, line)
			}
			panic(err)
		}
		EToV.Set(ind-1, 0, float64(n1-1))
		EToV.Set(ind-1, 1, float64(n2-1))
		EToV.Set(ind-1, 2, float64(n3-1))
	}
	return
}

func CalculateGeomAttributes(Nsd int, triIn3d bool) (NFaces int, bIs3D, bCoord3D, bElement3D, bTET bool) {
	bIs3D = (Nsd == 3)
	if bIs3D && !triIn3d {
		NFaces = 4 // Tetrahedra
		bCoord3D = true
		bElement3D = true
	} else {
		NFaces = 3 // Triangles
		bCoord3D = triIn3d
		bElement3D = false
	}
	// Triangles or Tetrahedra?
	bTET = bElement3D
	return
}

func getLine(reader *bufio.Reader) (line string) {
	var (
		err error
	)
	line, err = reader.ReadString('\n')
	if err != nil {
		if err == io.EOF {
			err = fmt.Errorf("early end of file")
		}
		panic(err)
	}
	line = line[:len(line)-1] // Strip away the newline
	return
}

func skipLines(n int, reader *bufio.Reader) {
	for i := 0; i < n; i++ {
		getLine(reader)
	}
}

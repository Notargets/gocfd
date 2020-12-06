package readfiles

import (
	"bufio"
	"fmt"
	"os"

	"github.com/notargets/gocfd/utils"
)

func ReadSU2(filename string, verbose bool) (K int, VX, VY utils.Vector, EToV, BCType utils.Matrix) {
	var (
		file   *os.File
		err    error
		reader *bufio.Reader
	)
	if verbose {
		fmt.Printf("Reading SU2 file named: %s\n", filename)
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
	BCType = ReadBCS(Nbcs, K, NFaces, reader)
	return
}

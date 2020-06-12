package DG2D

import (
	"bufio"
	"fmt"
	"io"
	"os"

	"github.com/notargets/gocfd/utils"
)

func ReadGambit2d(filename string) {
	var (
		file       *os.File
		err        error
		reader     *bufio.Reader
		VX, VY, VZ utils.Vector
	)
	fmt.Printf("Reading file named: %s\n", filename)
	if file, err = os.Open(filename); err != nil {
		panic(fmt.Errorf("unable to read file %s\n %s", filename, err))
	}
	defer file.Close()
	reader = bufio.NewReader(file)
	var line string
	// Skip firest six lines
	skipLines(6, reader)
	// Get dimensions
	/*
		Nv      // num nodes in mesh
		K       // num elements
		Nmats   // num material groups
		Nbcs    // num boundary groups
		Nsd;    // num space dimensions
	*/
	var Nv, K, Nmats, Nbcs, Nsd, dum, n int
	line = getLine(reader)
	if n, err = fmt.Sscanf(line, "%d %d %d %d %d %d", &Nv, &K, &Nmats, &Nbcs, &Nsd, &dum); err != nil || n < 5 {
		if err == nil && n < 5 {
			err = fmt.Errorf("read fewer than 5 dimensions, read %d, need 5\n, line: %s", n, line)
		}
		panic(err)
	}
	skipLines(2, reader)
	fmt.Printf("Nv = %d, K = %d\n", Nv, K)
	fmt.Printf("Nmats = %d, Nbcs = %d\n%d space dimensions\n", Nmats, Nbcs, Nsd)
	if Nsd > 3 || Nsd < 2 {
		panic("space dimensions not 2 or 3")
	}
	NFaces, bIs3D, bCoord3D, bElement3D, bTET := GetGeomAttributes(Nsd, false)
	fmt.Printf("NFaces = %d, bIs3D is %v, bCoord3D is %v, bElement3D is %v, bTET is %v\n",
		NFaces, bIs3D, bCoord3D, bElement3D, bTET)
	if bCoord3D {
		VX, VY, VZ = Read3DVertices(Nv, reader)
	} else {
		VX, VY = Read2DVertices(Nv, reader)
	}
	fmt.Println(VX.Print("VX"))
	fmt.Println(VY.Print("VY"))
	if bCoord3D {
		fmt.Println(VZ.Print("VZ"))
	}
}

func GetGeomAttributes(Nsd int, triIn3d bool) (NFaces int, bIs3D, bCoord3D, bElement3D, bTET bool) {
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

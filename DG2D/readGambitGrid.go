package DG2D

import (
	"bufio"
	"fmt"
	"io"
	"os"
)

func ReadGambit2d(filename string) {
	var (
		file   *os.File
		err    error
		reader *bufio.Reader
	)
	fmt.Printf("Reading file named: %s\n", filename)
	if file, err = os.Open(filename); err != nil {
		panic(fmt.Errorf("unable to read file %s\n %s", filename, err))
	}
	defer file.Close()
	reader = bufio.NewReader(file)
	var line string
	// Skip firest six lines
	for i := 0; i < 6; i++ {
		line = getLine(reader)
	}
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
	fmt.Printf("Nv = %d, K = %d\n", Nv, K)
	fmt.Printf("Nmats = %d, Nbcs = %d\n%d space dimensions\n", Nmats, Nbcs, Nsd)
	if Nsd > 3 || Nsd < 2 {
		panic("space dimensions not 2 or 3")
	}
	NFaces, bIs3D, bCoord3D, bElement3D, bTET := GetGeomAttributes(Nsd, false)
	fmt.Printf("NFaces = %d, bIs3D is %v, bCoord3D is %v, bElement3D is %v, bTET is %v\n",
		NFaces, bIs3D, bCoord3D, bElement3D, bTET)
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

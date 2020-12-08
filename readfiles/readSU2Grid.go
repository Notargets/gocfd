package readfiles

import (
	"bufio"
	"fmt"
	"os"
	"strings"

	"github.com/notargets/gocfd/types"

	"github.com/notargets/gocfd/utils"
)

// From here: https://su2code.github.io/docs_v7/Mesh-File/
type SU2ElementType uint8

const (
	ELType_LINE          SU2ElementType = 3
	ELType_Triangle                     = 5
	ELType_Quadrilateral                = 9
	ELType_Tetrahedral                  = 10
	ELType_Hexahedral                   = 12
	ELType_Prism                        = 13
	ELType_Pyramid                      = 14
)

func readBCs(reader *bufio.Reader) (BCEdges map[types.BCTAG][]types.EdgeInt) {
	var (
		nType   int
		v1, v2  int
		err     error
		prevInd int
	)
	NBCs := readNumber(reader)
	BCEdges = make(map[types.BCTAG][]types.EdgeInt, NBCs)
	for n := 0; n < NBCs; n++ {
		label := readLabel(reader)
		key := types.NewBCTAG(label)
		nEdges := readNumber(reader)
		if _, ok := BCEdges[key]; !ok {
			BCEdges[key] = make([]types.EdgeInt, nEdges)
		} else {
			// Add Nedges more storage to the map
			// This will end up appending duplicate tagged BCs to a common slice
			// For instance: Periodic BCs should come in pairs, so there should be Nedges x 2 edges
			prevInd = len(BCEdges[key])
			bsI := types.GrowSlice(BCEdges[key], prevInd+nEdges)
			BCEdges[key] = bsI.([]types.EdgeInt)
		}
		for i := 0; i < nEdges; i++ {
			line := getLine(reader)
			if _, err = fmt.Sscanf(line, "%d %d %d", &nType, &v1, &v2); err != nil {
				panic(err)
			}
			if SU2ElementType(nType) != ELType_LINE {
				panic("BCs should only contain line elements in 2D")
			}
			BCEdges[key][i+prevInd] = types.NewEdgeInt([2]int{v1, v2})
		}
	}
	return
}

func readVertices(reader *bufio.Reader) (VX, VY utils.Vector) {
	var (
		n    int
		x, y float64
		err  error
	)
	Nv := readNumber(reader)
	VX, VY = utils.NewVector(Nv), utils.NewVector(Nv)
	vxD, vyD := VX.Data(), VY.Data()
	for i := 0; i < Nv; i++ {
		line := getLine(reader)
		if n, err = fmt.Sscanf(line, "%f %f", &x, &y); err != nil {
			panic(err)
		}
		if n != 2 {
			panic("unable to read coordinates")
		}
		vxD[i], vyD[i] = x, y
	}
	return
}

func readElements(reader *bufio.Reader) (K int, EToV utils.Matrix) {
	var (
		n          int
		nType      int
		v1, v2, v3 int
		err        error
	)
	// EToV is K x 3
	K = readNumber(reader)
	EToV = utils.NewMatrix(K, 3)
	for k := 0; k < K; k++ {
		line := getLine(reader)
		if n, err = fmt.Sscanf(line, "%d %d %d %d", &nType, &v1, &v2, &v3); err != nil {
			panic(err)
		}
		if n != 4 {
			panic("unable to read vertices")
		}
		if SU2ElementType(nType) != ELType_Triangle {
			panic("unable to deal with non-triangular elements right now")
		}
		EToV.Set(k, 0, float64(v1))
		EToV.Set(k, 1, float64(v2))
		EToV.Set(k, 2, float64(v3))
	}
	return
}

func getToken(reader *bufio.Reader) (token string) {
	var (
		line string
		err  error
	)
	line = getLineNoComments(reader)
	ind := strings.Index(line, "=")
	if ind < 0 {
		err = fmt.Errorf("badly formed input line [%s], should have an =", line)
		panic(err)
	}
	token = line[ind+1:]
	return
}

func readLabel(reader *bufio.Reader) (label string) {
	var (
		err error
	)
	token := getToken(reader)
	if _, err = fmt.Sscanf(token, "%s", &label); err != nil {
		err = fmt.Errorf("unable to read label from token: [%s]", token)
		panic(err)
	}
	label = strings.Trim(label, " ")
	return
}

func readNumber(reader *bufio.Reader) (num int) {
	var (
		err error
	)
	token := getToken(reader)
	if _, err = fmt.Sscanf(token, "%d", &num); err != nil {
		err = fmt.Errorf("unable to read number from token: [%s]", token)
		panic(err)
	}
	return
}

func getLineNoComments(reader *bufio.Reader) (line string) {
	var ()
	for {
		line = strings.Trim(getLine(reader), " ")
		//fmt.Printf("line = [%s]\n", line)
		ind := strings.Index(line, "%")
		if ind < 0 || ind != 0 {
			return
		}
	}
}

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

	dimensionality := readNumber(reader)
	fmt.Printf("Read file with %d dimensional data...\n", dimensionality)
	K, EToV = readElements(reader)
	VX, VY = readVertices(reader)
	BCEdges := readBCs(reader)
	_ = BCEdges
	// TODO: Replace BCType with BCEdges above this call
	return
}

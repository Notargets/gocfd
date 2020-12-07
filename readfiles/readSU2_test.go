package readfiles

import (
	"bufio"
	"bytes"
	"fmt"
	"testing"

	"github.com/notargets/gocfd/types"

	"github.com/stretchr/testify/assert"
)

func TestReadSU2(t *testing.T) {
	{ // Test reading the file structure
		reader := bufio.NewReader(bytes.NewReader(inputFile))

		dim := readNumber(reader)
		assert.Equal(t, 2, dim)
		nelem := readNumber(reader)
		assert.Equal(t, 22, nelem)
		skipLines(22, reader)
		npts := readNumber(reader)
		assert.Equal(t, 18, npts)
		skipLines(18, reader)
		nmark := readNumber(reader)
		assert.Equal(t, 4, nmark)
		labels := []string{"periodic-left", "periodic-right", "top", "bottom"}
		nptsBC := []int{2, 2, 4, 4}
		for n := 0; n < nmark; n++ {
			mark := readLabel(reader)
			assert.Equal(t, labels[n], mark)
			nm := readNumber(reader)
			assert.Equal(t, nptsBC[n], nm)
			skipLines(nm, reader)
		}
	}
	{ // Test read elements and vertices
		reader := bufio.NewReader(bytes.NewReader(inputFile))
		_ = readNumber(reader)
		K, EToV := readElements(reader)
		assert.Equal(t, 22, K)
		assert.Equal(t, 17, int(EToV.At(K-1, 2)))
		VX, VY := readVertices(reader)
		Nv, _ := VX.Dims()
		assert.Equal(t, 18, Nv)
		Nv, _ = VY.Dims()
		assert.Equal(t, 18, Nv)
		assert.Equal(t, -7.100939331382065, VX.Data()[Nv-1])
		assert.Equal(t, 2.889910324036197, VY.Data()[Nv-1])
	}
	{ // Test read BCs
		reader := bufio.NewReader(bytes.NewReader(inputFile))
		_ = readNumber(reader)
		_, _ = readElements(reader)
		_, _ = readVertices(reader)
		BCEdges := readBCs(reader)
		_ = BCEdges
	}
}

func readBCs(reader *bufio.Reader) (BCEdges map[string][]types.EdgeNumber) {
	var (
		nType  int
		v1, v2 int
		err    error
	)
	NBCs := readNumber(reader)
	BCEdges = make(map[string][]types.EdgeNumber, NBCs)
	for n := 0; n < NBCs; n++ {
		label := readLabel(reader)
		if _, ok := BCEdges[label]; ok {
			err = fmt.Errorf("duplicate boundary condition found with label: [%s]", label)
			panic(err)
		}
		nEdges := readNumber(reader)
		BCEdges[label] = make([]types.EdgeNumber, nEdges)
		for i := 0; i < nEdges; i++ {
			line := getLine(reader)
			if _, err = fmt.Sscanf(line, "%d %d %d", &nType, &v1, &v2); err != nil {
				panic(err)
			}
			if SU2ElementType(nType) != ELType_LINE {
				panic("BCs should only contain line elements in 2D")
			}
			BCEdges[label][i] = types.NewEdgeNumber([2]int{v1, v2})
		}
	}
	return
}

var (
	inputFile = []byte(` %This is an example input file in SU2 format, output from gmsh
% Comments can appear outside of data areas
NDIME= 2
% Comments can appear outside of data areas
NELEM= 22
5 5 6 13 0
5 9 10 12 1
5 12 5 13 2
5 9 12 13 3
5 13 6 14 4
5 12 10 15 5
5 8 9 13 6
5 4 5 12 7
5 1 7 14 8
5 6 1 14 9
5 3 11 15 10
5 10 3 15 11
5 8 13 16 12
5 4 12 17 13
5 13 14 16 14
5 12 15 17 15
5 7 2 16 16
5 11 0 17 17
5 2 8 16 18
5 0 4 17 19
5 14 7 16 20
5 15 11 17 21
% Comments can appear outside of data areas
NPOIN= 18
-10 0 0
10 0 1
10 10 2
-10 10 3
-5.000000000004944 0 4
-1.231725832440134e-11 0 5
4.99999999999384 0 6
10 4.999999999992398 7
5.000000000004944 10 8
1.231725832440134e-11 10 9
-4.99999999999384 10 10
-10 5 11
-2.500000000008632 4.330127018915808 12
2.50000000000863 5.669872981084192 13
6.712741669205853 3.668411415814691 14
-6.712741669205681 6.331588584184096 15
7.100939331384343 7.110089675963254 16
-7.100939331382065 2.889910324036197 17
NMARK= 4
% Comments can appear outside of data areas
MARKER_TAG= periodic-left
% Comments can appear outside of data areas
MARKER_ELEMS= 2
3 3 11
3 11 0
% Comments can appear outside of data areas
MARKER_TAG= periodic-right
MARKER_ELEMS= 2
3 1 7
3 7 2
% Comments can appear outside of data areas
MARKER_TAG= top
MARKER_ELEMS= 4
3 2 8
3 8 9
3 9 10
3 10 3
MARKER_TAG= bottom
% Comments can appear outside of data areas
MARKER_ELEMS= 4
3 0 4
3 4 5
3 5 6
3 6 1
% Comments can appear outside of data areas
`)
)

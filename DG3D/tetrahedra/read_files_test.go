package tetrahedra

import (
	"testing"
)

func TestReadGambitNeutralFile(t *testing.T) {
	gf, err := ReadGambitNeutralFile("cube-partitioned.neu")
	if err != nil {
		panic(err)
	}
	gf.PrintMeshInfo()
}

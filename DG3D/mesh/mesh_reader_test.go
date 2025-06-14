package mesh

import (
	"testing"
)

func TestReadGambitNeutralFile(t *testing.T) {
	gf, err := ReadMeshFile("cube-partitioned.neu")
	if err != nil {
		panic(err)
	}
	gf.PrintStatistics()
}

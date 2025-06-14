package mesh

import (
	"github.com/stretchr/testify/assert"
	"testing"
)

func TestReadGambitNeutralFile(t *testing.T) {
	gf, err := ReadMeshFile("cube-partitioned.neu")
	if err != nil {
		panic(err)
	}
	gf.PrintStatistics()
	assert.Equal(t, 565, gf.NumElements)
	assert.Equal(t, 175, gf.NumVertices)
}

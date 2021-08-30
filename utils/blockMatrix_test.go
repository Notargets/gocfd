package utils

import (
	"fmt"
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestBlockMatrix(t *testing.T) {
	Bm := NewBlockMatrix(4)
	assert.True(t, Bm.M != nil)
	assert.Equal(t, 4, Bm.N)
	Bm.M[0][0] = NewMatrix(2, 2, []float64{
		0, 1,
		2, 3,
	})
	Bm.M[0][1] = NewMatrix(2, 2)
	Bm.M[1][1] = NewMatrix(2, 2)
	Bm.M[2][1] = NewMatrix(1, 1, []float64{100.})
	fmt.Printf(Bm.Print())
}

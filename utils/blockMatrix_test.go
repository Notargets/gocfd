package utils

import (
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestBlockMatrix(t *testing.T) {
	Bm := NewBlockMatrix(4, 2)
	assert.True(t, Bm.A != nil)
	assert.Equal(t, 4, Bm.N)
	assert.Equal(t, 2, Bm.NB)
	Bm.A[0][0] = NewMatrix(2, 2, []float64{
		0, 1,
		2, 3,
	})
	Bm.A[0][1] = NewMatrix(2, 2)
	Bm.A[1][1] = NewMatrix(2, 2)
	Bm.A[2][1] = NewMatrix(1, 1, []float64{100.})
	//fmt.Printf(Bm.Print())
}

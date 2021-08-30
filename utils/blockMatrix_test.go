package utils

import (
	"fmt"
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestBlockMatrix(t *testing.T) {
	Bm := NewBlockMatrix(4, 2)
	assert.True(t, Bm.M != nil)
	assert.Equal(t, 4, Bm.N)
	assert.Equal(t, 2, Bm.NB)
	fmt.Printf(Bm.Print())
}

package utils

import (
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestRanger(t *testing.T) {
	var (
		i1, i2 int
	)
	// Dimension parsing
	{
		i1, i2 = ParseDim(":", 10)
		assert.Equal(t, 0, i1)
		assert.Equal(t, 10, i2)
		i1, i2 = ParseDim(":5", 10)
		assert.Equal(t, 0, i1)
		assert.Equal(t, 5, i2)
		i1, i2 = ParseDim("5:5", 10)
		assert.Equal(t, 5, i1)
		assert.Equal(t, 6, i2)
		i1, i2 = ParseDim(4, 10)
		assert.Equal(t, 4, i1)
		assert.Equal(t, 5, i2)
		i1, i2 = ParseDim("2", 10)
		assert.Equal(t, 2, i1)
		assert.Equal(t, 3, i2)
	}
	// R2 indexing
	{
		my2d := NewR2(3, 4)
		index := my2d.Range(0, 0)
		assert.Equal(t, Index{0}, index)

		index = my2d.Range(0, ":")
		assert.Equal(t, Index{0, 1, 2, 3}, index)

		index = my2d.Range(0, ":")
		assert.Equal(t, Index{0, 1, 2, 3}, index)
		index = my2d.Range(1, ":")
		assert.Equal(t, Index{4, 5, 6, 7}, index)
		index = my2d.Range(":", 0)
		assert.Equal(t, Index{0, 4, 8}, index)
		index = my2d.Range(":", 1)
		assert.Equal(t, Index{1, 5, 9}, index)
	}
	// R3 indexing
	{
		my3d := NewR3(2, 3, 4)
		index := my3d.Range(0, 0, 0)
		assert.Equal(t, Index{0}, index)
		index = my3d.Range(0, ":", 0)
		assert.Equal(t, Index{0, 1, 2}, index)
		index = my3d.Range(":", 0, 0)
		assert.Equal(t, Index{0, 3}, index)
		index = my3d.Range(0, 0, ":")
		assert.Equal(t, Index{0, 6, 12, 18}, index)
		index = my3d.Range(1, 0, ":")
		assert.Equal(t, Index{3, 9, 15, 21}, index)
		index = my3d.Range(0, 1, ":")
		assert.Equal(t, Index{1, 7, 13, 19}, index)
	}
}

package DG2D

import (
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestRTElement_CalculateBasis(t *testing.T) {
	{
		rtb := NewRTBasis2DSimplex(3)
		assert.Equal(t, rtb.P, 3)
	}
}
package DG1D

import (
	"fmt"
	"testing"

	"github.com/stretchr/testify/assert"

	"gonum.org/v1/gonum/mat"

	"github.com/notargets/gocfd/utils"
)

func TestOperators(t *testing.T) {
	{
		var (
			A = utils.NewVector(4, []float64{3, 3, 1, 2})
			B = utils.NewVector(4, []float64{4, 4, 4, 4})
			C = utils.NewVector(4, []float64{5, -5, 5, 5})
		)
		D := Minmod(A, B, C)
		assert.Equal(t, utils.NewVector(4, []float64{3, 0, 1, 2}), D)
		fmt.Printf("D = \n%v\n", mat.Formatted(D, mat.Squeeze()))
	}
}

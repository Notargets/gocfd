package DG2D

import (
	"fmt"
	"strconv"
	"testing"
	"time"

	"github.com/stretchr/testify/assert"

	"github.com/notargets/gocfd/utils"

	"github.com/notargets/gocfd/DG1D"
)

func TestNewLagrangePolynomial1D(t *testing.T) {
	P := 2
	R := utils.NewVector(P+1, DG1D.LegendreZeros(P))

	lp1d := NewLagrangePolynomial1D(R, P, 0, 0)
	// Check the Lagrange property
	for j := 0; j < lp1d.Np; j++ {
		for i := 0; i < lp1d.Np; i++ {
			r := R.AtVec(i)
			psi := lp1d.getPolynomial(r, j)
			fmt.Printf("psi[%d,%d,%f]=%f\n", i, j, r, psi)
			if i == j {
				assert.InDeltaf(t, psi, 1, 0.000001, "")
			} else {
				assert.InDeltaf(t, psi, 0, 0.000001, "")
			}
		}
	}
}

func TestLagrangePolynomial2D(t *testing.T) {
	// Check the Lagrange Polynomial basis to verify the Lagrange property
	for P := 1; P < 7; P++ {
		R, S := NodesEpsilon(P)
		lp2d := NewLagrangePolynomialBasis2D(P, R, S)
		A := utils.NewMatrix(lp2d.Np, lp2d.Np)
		// Evaluate the j-th lagrange polynomial at (r,s)i
		// It should evaluate to 1 at each (r,s)i=j and 0 at (r,s)i!=j
		// In other words, this should be the diagonal matrix
		for j := 0; j < lp2d.Np; j++ {
			for i := 0; i < lp2d.Np; i++ {
				r, s := R.AtVec(i), S.AtVec(i)
				// fmt.Printf("psi[%d][%f,%f] = %f\n", j, r, s,
				// 	lp2d.GetPolynomialEvaluation(r, s, j))
				A.Set(i, j, lp2d.GetPolynomialEvaluation(r, s, j))
			}
		}
		checkIfUnitMatrix(t, A)
	}
}

func TestLagrangePolynomial(t *testing.T) {
	numSamples := 100
	rd := make([]float64, numSamples)
	xmin, xmax := -1., 1.
	fmin, fmax := -0.5, 1.25
	inc := (xmax - xmin) / float64(numSamples-1.)
	for i := 0; i < numSamples; i++ {
		rd[i] = xmin + float64(i)*inc
	}
	SamplesR := utils.NewVector(numSamples, rd)
	f := make([]float64, SamplesR.Len())
	var plot bool
	plot = false
	if plot {
		chart := utils.NewLineChart(1920, 1080, xmin, xmax, fmin, fmax)
		// TODO: Make a pluggable basis underneath the RT (and Lagrange) elements - Lagrange, Hesthaven, Spectral?
		var delay time.Duration
		Nmax := 4
		lineInc := 2. / float64(Nmax-2)
		lineColor := -1. // colormap goes from -1,1
		for n := 1; n < Nmax; n++ {
			Np := n + 1
			R := utils.NewVector(Np)
			inc = (xmax - xmin) / float64(Np-1)
			for i := 0; i < Np; i++ {
				R.DataP[i] = xmin + float64(i)*inc
			}
			lp := NewLagrangeBasis1D(R.DataP)
			_ = lp
			for j := 0; j < R.Len(); j++ {
				for i, r := range SamplesR.DataP {
					// f[i] = DG1D.JacobiP(utils.NewVector(1, []float64{r}), 0, 0, j)[0]
					f[i] = lp.BasisPolynomial([]float64{r}, j)[0]
				}
				if n == Nmax-1 && j == R.Len()-1 {
					delay = 120 * time.Second
				}
				name := "JacobiP[" + strconv.Itoa(n) + "," + strconv.Itoa(j) + "]"
				fmt.Printf("Chart Name: [%s], lineColor = %5.3f\n", name, lineColor)
				chart.Plot(delay, SamplesR.DataP, f, lineColor, name)
			}
			lineColor += lineInc
		}
	}

}

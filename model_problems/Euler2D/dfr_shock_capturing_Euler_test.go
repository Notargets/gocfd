package Euler2D

import (
	"fmt"
	"math"
	"testing"

	"github.com/stretchr/testify/assert"

	"github.com/notargets/gocfd/utils"

	"github.com/notargets/gocfd/DG2D"
)

func TestElementMean(t *testing.T) {
	var (
		N  = 2
		ip = ipDefault
	)
	ip.PolynomialOrder = N
	if !testing.Verbose() {
		return
	}
	// DFR := NewDFR2D(N, false, "../../DG2D/test_data/test_10tris_centered.
	// neu")
	meshFile := "../../DG2D/test_data/test_10tris_centered.neu"
	c := NewEuler(ip, meshFile, 1, false, false)
	// DG2D.SetTestFieldQ(c.DFR, DG2D.RADIAL1TEST, c.Q[0])
	// DG2D.SetTestFieldQ(c.DFR, DG2D.NORMALSHOCKTESTM12, c.Q[0])
	// DG2D.SetTestFieldQ(c.DFR, DG2D.NORMALSHOCKTESTM18, c.Q[0])
	// DG2D.SetTestFieldQ(c.DFR, DG2D.NORMALSHOCKTESTM2, c.Q[0])
	DG2D.SetTestFieldQ(c.DFR, DG2D.NORMALSHOCKTESTM5, c.Q[0])
	// DG2D.SetTestFieldQ(c.DFR, DG2D.FIXEDVORTEXTEST, c.Q[0])
	// DG2D.SetTestFieldQ(c.DFR, DG2D.INTEGERTEST, c.Q[0])
	_, KMax := c.Q[0][0].Dims()
	var QMean [4]utils.Vector
	for n := 0; n < 4; n++ {
		QMean[n] = utils.NewVector(KMax)
	}
	c.UpdateElementMean(c.Q[0], QMean)
	// for n := 0; n < 4; n++ {
	// 	fmt.Printf("QMean[%d]=", n)
	// 	c.RK.QMean[0][n].Transpose().Print("")
	// }
	assert.InDeltaSlicef(t,
		[]float64{1, 1, 1.40545, 4.28205, 4.6875, 1, 1, 1.40545, 4.28205, 4.6875},
		QMean[0].DataP, 1e-4, "")
	assert.InDeltaSlicef(t,
		[]float64{5, 5, 5, 5, 5, 5, 5, 5, 5, 5}, QMean[1].DataP, 1e-4, "")
	assert.InDeltaSlicef(t,
		[]float64{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, QMean[2].DataP, 1e-4, "")
	assert.InDeltaSlicef(t,
		[]float64{15, 15, 19.32477, 50.00856, 54.33333, 15, 15, 19.32477,
			50.00856, 54.33333},
		QMean[3].DataP, 1e-4, "")
}

func TestPlotShockTemperedInterpolation(t *testing.T) {
	var (
		N               = 2
		shockedElements = []int{2, 3, 7, 8}
		ip              = ipDefault
	)
	ip.PolynomialOrder = N
	if !testing.Verbose() {
		return
	}
	// DFR := NewDFR2D(N, false, "../../DG2D/test_data/test_10tris_centered.
	// neu")
	meshFile := "../../DG2D/test_data/test_10tris_centered.neu"
	c := NewEuler(ip, meshFile, 1, false, false)
	// DG2D.SetTestFieldQ(c.DFR, DG2D.RADIAL1TEST, c.Q[0])
	// DG2D.SetTestFieldQ(c.DFR, DG2D.NORMALSHOCKTESTM12, c.Q[0])
	// DG2D.SetTestFieldQ(c.DFR, DG2D.NORMALSHOCKTESTM18, c.Q[0])
	// DG2D.SetTestFieldQ(c.DFR, DG2D.NORMALSHOCKTESTM2, c.Q[0])
	DG2D.SetTestFieldQ(c.DFR, DG2D.NORMALSHOCKTESTM5, c.Q[0])
	// DG2D.SetTestFieldQ(c.DFR, DG2D.FIXEDVORTEXTEST, c.Q[0])
	// DG2D.SetTestFieldQ(c.DFR, DG2D.INTEGERTEST, c.Q[0])
	var QInterp [4]utils.Matrix
	for n := 0; n < 4; n++ {
		QInterp[n] = c.DFR.GraphInterp.Mul(c.Q[0][n])
	}
	c.modulateQInterp(c.Q[0], QInterp, c.ShockFinder)
	for _, k := range shockedElements {
		getElementInterpolationStats(c, QInterp, c.Q[0], k)
	}
	DG2D.PlotField(QInterp[0].Transpose().DataP, c.DFR.GraphMesh, 0.0, 0.0)
}

func getElementInterpolationStats(c *Euler, QInterp, Q [4]utils.Matrix, k int) {
	NpInt, _ := QInterp[0].Dims()
	var UDelta, UMean, UMin, UMax [4]float64
	UMean = c.getElementMean(Q, k)
	UMin, UMax = getElementMinMax(Q, k)
	fmt.Printf("UMean[%d] = {%.1f, %.1f, %.1f, %.1f}, ",
		k, UMean[0], UMean[1], UMean[2], UMean[3])
	fmt.Printf("UMin[%d] = {%.1f, %.1f, %.1f, %.1f}, ",
		k, UMin[0], UMin[1], UMin[2], UMin[3])
	fmt.Printf("UMax[%d] = {%.1f, %.1f, %.1f, %.1f} ",
		k, UMax[0], UMax[1], UMax[2], UMax[3])
	fmt.Println()
	n := 0
	fmt.Printf("Interp UDelta[%d]: ", k)
	for i := 0; i < NpInt; i++ {
		val := QInterp[n].At(i, k)
		UDelta[n] = 100. * (val - UMean[n]) / (UMean[n] + 1e-10)
		fmt.Printf("%.1f%%,", UDelta[n])
	}
	fmt.Println()
	fmt.Printf("Interp OverMax[%d]: ", k)
	for i := 0; i < NpInt; i++ {
		val := QInterp[n].At(i, k)
		delta := val - UMax[n]
		if delta > 0.0 {
			fmt.Printf("%.1f%%,", 100.*delta/UMean[n])
		}
	}
	fmt.Println()
	fmt.Printf("Interp UnderMin[%d]: ", k)
	for i := 0; i < NpInt; i++ {
		val := QInterp[n].At(i, k)
		delta := UMin[n] - val
		if delta > 0.0 {
			fmt.Printf("%.1f%%,", 100.*delta/UMean[n])
		}
	}
	fmt.Println()
	fmt.Println()
}

func getElementMinMax(Q [4]utils.Matrix, k int) (Min, Max [4]float64) {
	var (
		Np, _ = Q[0].Dims()
	)
	for n := 0; n < 4; n++ {
		Min[n] = math.MaxFloat64
		Max[n] = -math.MaxFloat64
		for i := 0; i < Np; i++ {
			Min[n] = min(Min[n], Q[n].At(i, k))
			Max[n] = max(Max[n], Q[n].At(i, k))
		}
	}
	return
}

func TestPlotEntropyInterpolation(t *testing.T) {
	var (
		N  = 2
		ip = ipDefault
	)
	ip.PolynomialOrder = N
	if !testing.Verbose() {
		return
	}
	meshFile := "../../DG2D/test_data/test_10tris_centered.neu"
	c := NewEuler(ip, meshFile, 1, false, false)
	// DG2D.SetTestFieldQ(c.DFR, DG2D.RADIAL1TEST, c.Q[0])
	// DG2D.SetTestFieldQ(c.DFR, DG2D.NORMALSHOCKTESTM12, c.Q[0])
	// DG2D.SetTestFieldQ(c.DFR, DG2D.NORMALSHOCKTESTM18, c.Q[0])
	// DG2D.SetTestFieldQ(c.DFR, DG2D.NORMALSHOCKTESTM2, c.Q[0])
	DG2D.SetTestFieldQ(c.DFR, DG2D.NORMALSHOCKTESTM5, c.Q[0])
	// DG2D.SetTestFieldQ(c.DFR, DG2D.FIXEDVORTEXTEST, c.Q[0])
	// DG2D.SetTestFieldQ(c.DFR, DG2D.INTEGERTEST, c.Q[0])
	// fMGraph := c.DFR.GraphInterp.Mul(c.Q[0][3]).Transpose()
	SwitchToEntropyVariables(c.Q[0], c.FSFar.Gamma)
	var QInterp [4]utils.Matrix
	for n := 0; n < 4; n++ {
		QInterp[n] = c.DFR.GraphInterp.Mul(c.Q[0][n])
	}
	SwitchToConservedVariables(QInterp, c.FSFar.Gamma)
	DG2D.PlotField(QInterp[0].Transpose().DataP, c.DFR.GraphMesh, 0.0, 0.0)
}

package DG2D

import (
	"fmt"
	"testing"
)

func TestInterpolationVariousFields(t *testing.T) {
	var (
		NMin = 7
		NMax = 7
	)
	for N := NMin; N <= NMax; N++ {
		fmt.Printf("ORDER: %d Element Test\n-----------------------\n", N)
		// angle := 210.
		angle := 82.8
		dfr := CreateEquiTriMesh(N, angle)
		// for _, tf := range []TestField{NORMALSHOCKTESTM12, NORMALSHOCKTESTM2,
		// 	NORMALSHOCKTESTM5, FIXEDVORTEXTEST, RADIAL1TEST, RADIAL2TEST,
		// 	RADIAL3TEST, RADIAL4TEST} {
		for _, tf := range []TestField{FIXEDVORTEXTEST, NORMALSHOCKTESTM5} {
			Np := dfr.SolutionElement.Np
			X, Y := dfr.SolutionX.DataP, dfr.SolutionY.DataP
			QSol := QFromField(setTestField(X, Y, tf), Np)
			fmt.Printf("%s Interpolation\n", tf.String())
			stats := GetInterpolationAccuracy(dfr, QSol, dfr.FluxEdgeInterp, tf)
			printRMSError(stats)
			printRMSError(stats)
			fmt.Printf("---------------------\n\n")
			// Coeffs := GetRTCoefficients(dfr, tf)
			// QSol.Print("QSol Orig")
			// for _, iterCount := range []int{10, 100, 1000} {
			// 	QSolMod := ModulateInternalField(dfr, QSol, Nu, p, iterCount)
			// 	QSolMod.Print("QSolMod." + strconv.Itoa(iterCount))
			// }
		}
	}
}

package DG2D

import (
	"math"
	"testing"

	"github.com/notargets/gocfd/utils"
)

func TestCorrectModalTransferWithActualBasis(t *testing.T) {
	N := 3
	t.Logf("=== CORRECT MODAL TRANSFER WITH ACTUAL 2D BASIS ===")

	// Create elements
	le := NewLagrangeElement2D(N)
	rt := NewRTElement(N+1, SimplexRTBasis, GaussEdgePoints)

	t.Logf("DG element: %d interior points", le.Np)

	// Test polynomial: f(r,s) = 1 + r + s
	testPoly := func(r, s float64) float64 {
		return 1.0 + r + s
	}

	// Evaluate at DG interior points
	interiorValues := make([]float64, le.Np)
	for i := 0; i < le.Np; i++ {
		r, s := le.R.AtVec(i), le.S.AtVec(i)
		interiorValues[i] = testPoly(r, s)
	}

	// Convert to modal coefficients using existing system
	interiorMatrix := utils.NewMatrix(len(interiorValues), 1, interiorValues)
	modalCoeffs2D := le.JB2D.Vinv.Mul(interiorMatrix)

	t.Logf("2D modal coefficients from your system:")
	for i := 0; i < modalCoeffs2D.Rows(); i++ {
		order := le.JB2D.Order2DAtJ[i]
		coeff := modalCoeffs2D.At(i, 0)
		t.Logf("  Mode (%d,%d): %.6f", order[0], order[1], coeff)
	}

	// Verify modal coefficients by reconstruction
	t.Logf("\nVerifying modal coefficients by reconstruction:")
	testPoints := [][]float64{
		{-0.5, -0.5},
		{0.0, 0.0},
		{0.3, 0.2},
	}

	for _, pt := range testPoints {
		r, s := pt[0], pt[1]
		exact := testPoly(r, s)

		// Reconstruct using the modal coefficients and your 2D basis
		reconstructed := reconstructUsing2DBasis(le, modalCoeffs2D, r, s)

		error := math.Abs(exact - reconstructed)
		t.Logf("  Point (%.1f,%.1f): exact=%.6f, reconstructed=%.6f, error=%.2e",
			r, s, exact, reconstructed, error)
	}

	// Test modal transfer for each edge
	for edgeNum := 0; edgeNum < 3; edgeNum++ {
		t.Logf("\n--- EDGE %d MODAL TRANSFER ---", edgeNum)

		// Build correct modal transfer matrix using actual understanding of 2D basis
		transferMatrix := buildActualModalTransferMatrix(le, rt, edgeNum)

		// Transfer to 1D modes
		modalCoeffs1D := transferMatrix.Mul(modalCoeffs2D)

		t.Logf("1D modal coefficients:")
		for i := 0; i < modalCoeffs1D.Rows(); i++ {
			coeff := modalCoeffs1D.At(i, 0)
			t.Logf("  1D Mode %d: %.6f", i, coeff)
		}

		// Convert to nodal values using your 1D Jacobi basis
		edgePoints := getEdgePoints(rt, edgeNum)
		edge1DCoords := make([]float64, len(edgePoints))
		for i, pt := range edgePoints {
			edge1DCoords[i] = mapToEdge1D(pt.R, pt.S, edgeNum)
		}

		edgeOrder := rt.NpEdge - 1
		edge1DCoordVec := utils.NewVector(len(edge1DCoords), edge1DCoords)
		jb1d := NewJacobiBasis1D(edgeOrder, edge1DCoordVec)
		Vedge1D := jb1d.Vandermonde1D()

		edgeMatrix := Vedge1D.Mul(modalCoeffs1D)

		// Compare with exact values
		maxError := 0.0
		for i, pt := range edgePoints {
			exact := testPoly(pt.R, pt.S)
			computed := edgeMatrix.At(i, 0)
			error := math.Abs(exact - computed)
			if error > maxError {
				maxError = error
			}
			t.Logf("  Point %d: exact=%.6f, computed=%.6f, error=%.2e",
				i, exact, computed, error)
		}

		t.Logf("Maximum error: %.2e", maxError)

		if maxError < 1e-10 {
			t.Logf("✓ Modal transfer SUCCESS for edge %d", edgeNum)
		} else {
			t.Logf("✗ Modal transfer FAILED for edge %d", edgeNum)
		}
	}
}

func reconstructUsing2DBasis(le *LagrangeElement2D, modalCoeffs utils.Matrix, r, s float64) float64 {
	// Reconstruct the polynomial at point (r,s) using modal coefficients and 2D basis
	result := 0.0

	RVec := utils.NewVector(1, []float64{r})
	SVec := utils.NewVector(1, []float64{s})

	for i := 0; i < modalCoeffs.Rows(); i++ {
		order := le.JB2D.Order2DAtJ[i]
		ii, jj := order[0], order[1]
		coeff := modalCoeffs.At(i, 0)

		// Evaluate 2D basis function at (r,s)
		basisValues := le.JB2D.Simplex2DP(RVec, SVec, ii, jj)
		basisValue := basisValues[0]

		result += coeff * basisValue
	}

	return result
}

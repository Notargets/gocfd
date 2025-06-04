package DG3D

import (
	"math"
	"testing"
)

func TestWilliamsShunnJamesonCubature_WeightScalingDiagnostic(t *testing.T) {
	// This test helps diagnose the weight scaling issue

	t.Log("=== Weight Scaling Diagnostic ===")
	t.Log("Reference tet vertices: (-1,-1,-1), (1,-1,-1), (-1,1,-1), (-1,-1,1)")
	t.Log("Reference tet volume: 4/3")
	t.Log("")

	// For each order, check the raw sum before and after the 4/3 scaling
	for order := 1; order <= 5; order++ {
		_, _, _, w := WilliamsShunnJamesonCubature(order)

		sum := 0.0
		for _, weight := range w {
			sum += weight
		}

		// Calculate what the unscaled sum would be
		unscaledSum := sum / (4.0 / 3.0)

		t.Logf("Order %d:", order)
		t.Logf("  Number of points: %d", len(w))
		t.Logf("  Sum of weights: %.15f", sum)
		t.Logf("  Unscaled sum (divided by 4/3): %.15f", unscaledSum)

		// Check if the unscaled sum is close to standard values
		if math.Abs(unscaledSum-1.0) < 1e-10 {
			t.Logf("  -> Weights appear to be for unit simplex (sum to 1) then scaled by 4/3")
		} else if math.Abs(unscaledSum-0.75) < 1e-10 {
			t.Logf("  -> Weights appear to be for a different scaling")
		}

		// For orders with issues, calculate the correction factor
		if math.Abs(sum-4.0/3.0) > 1e-14 {
			correctionFactor := (4.0 / 3.0) / sum
			t.Logf("  -> To fix: multiply all weights by %.15f", correctionFactor)
		}

		t.Log("")
	}

	// Check the specific published values for order 3
	t.Log("=== Order 3 Published Values Check ===")
	_, _, _, w3 := WilliamsShunnJamesonCubature(3)

	// The published weights from Williams-Shunn-Jameson for order 3
	wv_published := (25.0 - 5.0*math.Sqrt(5.0)) / 480.0
	we_published := (25.0 + 5.0*math.Sqrt(5.0)) / 480.0

	t.Logf("Published vertex weight: %.15f", wv_published)
	t.Logf("Published edge weight: %.15f", we_published)
	t.Logf("Sum of published weights: %.15f", 4*wv_published+6*we_published)

	// Check what we actually get
	t.Logf("Actual vertex weight: %.15f", w3[0]/(4.0/3.0))
	t.Logf("Actual edge weight: %.15f", w3[4]/(4.0/3.0))

	// The issue appears to be that the published weights sum to a different value
	publishedSum := 4*wv_published + 6*we_published
	t.Logf("Published weights sum to: %.15f (not 1!)", publishedSum)
	t.Logf("This explains why scaling by 4/3 gives: %.15f", publishedSum*4.0/3.0)

	// Additional analysis
	t.Log("")
	t.Log("=== Summary ===")
	t.Log("The Williams-Shunn-Jameson published weights are for a specific reference element.")
	t.Log("For orders 1-2: weights sum to 1 and are correctly scaled by 4/3")
	t.Log("For orders 3-5: weights do NOT sum to 1 before scaling, causing incorrect totals")
	t.Log("The implementation appears to incorrectly apply the 4/3 scaling to weights that")
	t.Log("were not designed to sum to 1.")
}

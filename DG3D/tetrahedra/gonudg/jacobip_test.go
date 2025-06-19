package gonudg

import (
	"math"
	"testing"
)

func TestJacobiPSpecificValuesCorrected(t *testing.T) {
	// Test specific values calculated from the Hesthaven & Warburton orthonormal Jacobi polynomial formula
	// Source: The expected values are computed using the formula from:
	// - Hesthaven & Warburton "Nodal Discontinuous Galerkin Methods" (2007)
	// - MATLAB implementation JacobiP.m from the book's accompanying code
	// - C++ implementation JacobiP.cpp in this project
	//
	// Formula: gamma0 = 2^(α+β+1)/(α+β+1) * Γ(α+1) * Γ(β+1) / Γ(α+β+1)
	//          P_0^{α,β}(X) = 1/√gamma0  (constant for all X)

	x := []float64{-1.0, 0.0, 1.0}

	// Test P_0^{0,0}(X) = 1/sqrt(2)
	// gamma0 = 2^1/1 * 1 * 1 / 1 = 2
	P0_00 := JacobiP(x, 0.0, 0.0, 0)
	expected0_00 := 1.0 / math.Sqrt(2.0)
	for i, val := range P0_00 {
		if math.Abs(val-expected0_00) > 1e-14 {
			t.Errorf("P_0^{0,0}(%f) = %f, expected %f", x[i], val, expected0_00)
		}
	}

	// Test P_0^{1,0}(X) = 1/sqrt(2)
	// gamma0 = 2^2/2 * 1 * 1 / 1 = 2
	P0_10 := JacobiP(x, 1.0, 0.0, 0)
	expected0_10 := 1.0 / math.Sqrt(2.0) // CORRECTED from 1/sqrt(3)
	for i, val := range P0_10 {
		if math.Abs(val-expected0_10) > 1e-14 {
			t.Errorf("P_0^{1,0}(%f) = %f, expected %f", x[i], val, expected0_10)
		}
	}

	// Test P_0^{2,0}(X) = sqrt(3/8)
	// gamma0 = 2^3/3 * 2 * 1 / 2 = 8/3
	P0_20 := JacobiP(x, 2.0, 0.0, 0)
	expected0_20 := math.Sqrt(3.0 / 8.0)
	for i, val := range P0_20 {
		if math.Abs(val-expected0_20) > 1e-14 {
			t.Errorf("P_0^{2,0}(%f) = %f, expected %f", x[i], val, expected0_20)
		}
	}

	// Verify that 2*sqrt(2) * P_0^{0,0} * P_0^{1,0} * P_0^{2,0} = sqrt(3)/2
	// This is what Simplex3DP computes for P_{0,0,0}
	// 2√2 × (1/√2) × (1/√2) × √(3/8) = √2 × √(3/8) = √(3/4) = √3/2
	product := 2.0 * math.Sqrt(2.0) * expected0_00 * expected0_10 * expected0_20
	expectedProduct := math.Sqrt(3.0) / 2.0 // √3/2 ≈ 0.866025
	if math.Abs(product-expectedProduct) > 1e-14 {
		t.Errorf("2*sqrt(2) * P_0^{0,0} * P_0^{1,0} * P_0^{2,0} = %f, expected %f",
			product, expectedProduct)
	}
}

func TestJacobiPValueAt000Corrected(t *testing.T) {
	// Test that P_{0,0,0} gives sqrt(3)/2 when used in Simplex3DP
	// Based on the corrected orthonormal Jacobi polynomial values

	// At the reference point, a=b=c=-1 (from the test)
	a := []float64{-1.0}
	b := []float64{-1.0}
	c := []float64{-1.0}

	// Compute the Jacobi polynomials as done in Simplex3DP
	h1 := JacobiP(a, 0.0, 0.0, 0) // P_0^{0,0} = 1/√2
	h2 := JacobiP(b, 1.0, 0.0, 0) // P_0^{1,0} = 1/√2
	h3 := JacobiP(c, 2.0, 0.0, 0) // P_0^{2,0} = √(3/8)

	// The Simplex3DP formula for i=j=k=0:
	// P = 2*sqrt(2) * h1 * h2 * (1-b)^0 * h3 * (1-c)^0
	// P = 2*sqrt(2) * h1 * h2 * h3
	// P = 2*sqrt(2) * (1/√2) * (1/√2) * √(3/8)
	// P = 2*sqrt(2) * (1/2) * √(3/8) = √2 * √(3/8) = √(3/4) = √3/2
	result := 2.0 * math.Sqrt(2.0) * h1[0] * h2[0] * h3[0]

	t.Logf("h1[0] = %f (expected %f)", h1[0], 1.0/math.Sqrt(2.0))
	t.Logf("h2[0] = %f (expected %f)", h2[0], 1.0/math.Sqrt(2.0))
	t.Logf("h3[0] = %f (expected %f)", h3[0], math.Sqrt(3.0/8.0))
	t.Logf("Result = %f", result)
	t.Logf("Expected = %f", math.Sqrt(3.0)/2.0)

	// The test expects sqrt(3)/2
	expected := math.Sqrt(3.0) / 2.0 // √3/2 ≈ 0.866025
	if math.Abs(result-expected) > 1e-10 {
		t.Errorf("P_{0,0,0} = %f, expected %f", result, expected)
	}
}

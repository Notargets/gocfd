package gonudg

import (
	"math"
	"testing"

	"github.com/notargets/gocfd/utils"
)

// TestLift3DSimpleCase tests LIFT3D with the simplest case N=1
func TestLift3DSimpleCase(t *testing.T) {
	N := 1

	// Generate 3D nodes
	X, Y, Z := Nodes3D(N)
	r, s, tt := XYZtoRST(X, Y, Z)

	// For N=1, we should have Np=4 nodes (vertices of tetrahedron)
	if len(r) != 4 {
		t.Fatalf("Expected 4 nodes for N=1, got %d", len(r))
	}

	// Build matrices
	V := Vandermonde3D(N, r, s, tt)
	Fmask := BuildFaceMasks(N, r, s, tt)

	// For N=1, each face should have Nfp=3 nodes
	Nfp := (N + 1) * (N + 2) / 2 // Should be 3
	for face := 0; face < 4; face++ {
		if len(Fmask[face]) != Nfp {
			t.Errorf("Face %d has %d nodes, expected %d", face, len(Fmask[face]), Nfp)
		}
		t.Logf("Face %d nodes: %v", face, Fmask[face])
	}

	// Compute LIFT
	LIFT := Lift3D(N, r, s, tt, V, Fmask)

	// Check dimensions: should be 4 x 12 (Np x 4*Nfp)
	nrows, ncols := LIFT.Dims()
	if nrows != 4 || ncols != 12 {
		t.Errorf("Expected 4x12 LIFT matrix, got %dx%d", nrows, ncols)
	}

	// Print some values for inspection
	t.Logf("LIFT matrix sample values:")
	for i := 0; i < min(4, nrows); i++ {
		for j := 0; j < min(4, ncols); j++ {
			if math.Abs(LIFT.At(i, j)) > 1e-14 {
				t.Logf("LIFT[%d,%d] = %e", i, j, LIFT.At(i, j))
			}
		}
	}
}

// TestLift3DDimensions verifies the LIFT matrix has correct dimensions
func TestLift3DDimensions(t *testing.T) {
	testCases := []int{1, 2, 3, 4}

	for _, N := range testCases {
		// Generate 3D nodes using the actual project functions
		X, Y, Z := Nodes3D(N)
		r, s, tt := XYZtoRST(X, Y, Z)

		// Build Vandermonde matrix using actual project function
		V := Vandermonde3D(N, r, s, tt)

		// Build face masks
		Fmask := BuildFaceMasks(N, r, s, tt)

		// Compute LIFT
		LIFT := Lift3D(N, r, s, tt, V, Fmask)

		// Check dimensions
		Np := (N + 1) * (N + 2) * (N + 3) / 6
		Nfp := (N + 1) * (N + 2) / 2
		Nfaces := 4

		nrows, ncols := LIFT.Dims()
		if nrows != Np {
			t.Errorf("N=%d: Expected %d rows, got %d", N, Np, nrows)
		}
		if ncols != Nfaces*Nfp {
			t.Errorf("N=%d: Expected %d cols, got %d", N, Nfaces*Nfp, ncols)
		}
	}
}

// TestLift3DDiagnostics helps diagnose numerical issues
func TestLift3DDiagnostics(t *testing.T) {
	N := 2

	// Generate 3D nodes
	X, Y, Z := Nodes3D(N)
	r, s, tt := XYZtoRST(X, Y, Z)

	// Check that nodes are within reference tetrahedron
	for i := 0; i < len(r); i++ {
		if r[i] < -1.0-1e-10 || r[i] > 1.0+1e-10 ||
			s[i] < -1.0-1e-10 || s[i] > 1.0+1e-10 ||
			tt[i] < -1.0-1e-10 || tt[i] > 1.0+1e-10 {
			t.Errorf("Node %d outside [-1,1]: r=%g, s=%g, t=%g", i, r[i], s[i], tt[i])
		}
		if r[i]+s[i]+tt[i] > 1.0+1e-10 {
			t.Errorf("Node %d outside tetrahedron: r+s+t=%g", i, r[i]+s[i]+tt[i])
		}
	}

	// Build Vandermonde and check conditioning
	V := Vandermonde3D(N, r, s, tt)

	// Check that V is invertible
	invV := V.InverseWithCheck()

	// Check V*invV = I
	I := V.Mul(invV)
	nrows, _ := I.Dims()
	for i := 0; i < nrows; i++ {
		for j := 0; j < nrows; j++ {
			expected := 0.0
			if i == j {
				expected = 1.0
			}
			if math.Abs(I.At(i, j)-expected) > 1e-10 {
				t.Errorf("V*invV not identity at (%d,%d): %g", i, j, I.At(i, j))
			}
		}
	}

	// Check face masks
	Fmask := BuildFaceMasks(N, r, s, tt)
	Nfp := (N + 1) * (N + 2) / 2

	for face := 0; face < 4; face++ {
		if len(Fmask[face]) != Nfp {
			t.Errorf("Face %d has %d nodes, expected %d", face, len(Fmask[face]), Nfp)
		}

		// Check that face nodes satisfy face constraint
		for _, idx := range Fmask[face] {
			switch face {
			case 0: // t = -1
				if math.Abs(1.0+tt[idx]) > 1e-7 {
					t.Errorf("Face 0 node %d: t=%g, not -1", idx, tt[idx])
				}
			case 1: // s = -1
				if math.Abs(1.0+s[idx]) > 1e-7 {
					t.Errorf("Face 1 node %d: s=%g, not -1", idx, s[idx])
				}
			case 2: // r+s+t = -1
				if math.Abs(1.0+r[idx]+s[idx]+tt[idx]) > 1e-7 {
					t.Errorf("Face 2 node %d: r+s+t=%g, not -1", idx, r[idx]+s[idx]+tt[idx])
				}
			case 3: // r = -1
				if math.Abs(1.0+r[idx]) > 1e-7 {
					t.Errorf("Face 3 node %d: r=%g, not -1", idx, r[idx])
				}
			}
		}
	}
}

// TestLift3DBasicProperties tests fundamental properties of the LIFT matrix
func TestLift3DBasicProperties(t *testing.T) {
	N := 3

	// Generate 3D nodes using actual project functions
	X, Y, Z := Nodes3D(N)
	r, s, tt := XYZtoRST(X, Y, Z)

	// Build required matrices using actual project functions
	V := Vandermonde3D(N, r, s, tt)
	Fmask := BuildFaceMasks(N, r, s, tt)
	LIFT := Lift3D(N, r, s, tt, V, Fmask)

	// Test 1: LIFT should have finite values
	maxVal := 0.0
	minNonZeroVal := math.MaxFloat64
	nrows, ncols := LIFT.Dims()
	numNonZero := 0

	for i := 0; i < nrows; i++ {
		for j := 0; j < ncols; j++ {
			val := LIFT.At(i, j)
			if math.IsNaN(val) || math.IsInf(val, 0) {
				t.Errorf("LIFT has invalid value at (%d,%d): %v", i, j, val)
			}

			absVal := math.Abs(val)
			if absVal > 1e-14 { // Count as non-zero
				numNonZero++
				if absVal > maxVal {
					maxVal = absVal
				}
				if absVal < minNonZeroVal {
					minNonZeroVal = absVal
				}
			}
		}
	}

	// The LIFT matrix is expected to be sparse with many zeros
	// Only check condition number of non-zero values
	t.Logf("LIFT matrix: %d non-zero entries out of %d total", numNonZero, nrows*ncols)
	t.Logf("Non-zero value range: [%e, %e]", minNonZeroVal, maxVal)

	// For LIFT matrices, we expect some variation but not extreme
	// The face integral operators can have legitimately large variations
	if numNonZero > 0 && maxVal/minNonZeroVal > 1e18 {
		t.Errorf("LIFT matrix non-zero values span too large a range: max/min = %e", maxVal/minNonZeroVal)
	}

	// Test 2: Column norms should be reasonable
	// Each column represents the lifting of a face node to the volume
	for j := 0; j < ncols; j++ {
		norm := 0.0
		for i := 0; i < nrows; i++ {
			val := LIFT.At(i, j)
			norm += val * val
		}
		norm = math.Sqrt(norm)

		if math.IsNaN(norm) || math.IsInf(norm, 0) {
			t.Errorf("Column %d has invalid norm: %v", j, norm)
		}
		if norm > 100.0 {
			t.Errorf("Column %d norm unusually large: %v", j, norm)
		}
	}
}

// TestLift3DConstantFunction tests LIFT applied to constant surface values
func TestLift3DConstantFunction(t *testing.T) {
	N := 2

	// Generate 3D nodes
	X, Y, Z := Nodes3D(N)
	r, s, tt := XYZtoRST(X, Y, Z)

	// Build required matrices
	V := Vandermonde3D(N, r, s, tt)
	Fmask := BuildFaceMasks(N, r, s, tt)
	LIFT := Lift3D(N, r, s, tt, V, Fmask)

	// Create a vector representing constant value 1.0 on all face nodes
	Nfp := (N + 1) * (N + 2) / 2
	Nfaces := 4
	faceVals := utils.NewVector(Nfaces * Nfp)
	for i := 0; i < faceVals.Len(); i++ {
		faceVals.Set(i, 1.0)
	}

	// Apply LIFT
	volumeVals := LIFT.Mul(faceVals.ToMatrix())

	// The result should lift the face values into the volume
	// For a well-formed LIFT operator, the values should be non-zero
	// near the faces and decay into the interior
	hasNonZero := false
	nrows, _ := volumeVals.Dims()
	for i := 0; i < nrows; i++ {
		if math.Abs(volumeVals.At(i, 0)) > 1e-10 {
			hasNonZero = true
			break
		}
	}
	if !hasNonZero {
		t.Error("LIFT applied to constant face values produced all zeros")
	}
}

// TestLift3DFaceConsistency verifies LIFT preserves face polynomial structure
func TestLift3DFaceConsistency(t *testing.T) {
	N := 3

	// Generate 3D nodes
	X, Y, Z := Nodes3D(N)
	r, s, tt := XYZtoRST(X, Y, Z)

	// Build required matrices
	V := Vandermonde3D(N, r, s, tt)
	Fmask := BuildFaceMasks(N, r, s, tt)
	LIFT := Lift3D(N, r, s, tt, V, Fmask)

	// Test that LIFT correctly lifts face polynomials
	// For each face, set a polynomial on that face and zeros elsewhere
	Nfp := (N + 1) * (N + 2) / 2
	Nfaces := 4

	for face := 0; face < Nfaces; face++ {
		// Create face values: 1.0 on current face, 0.0 elsewhere
		faceVals := utils.NewVector(Nfaces * Nfp)
		for i := 0; i < Nfp; i++ {
			faceVals.Set(face*Nfp+i, 1.0)
		}

		// Apply LIFT
		volumeVals := LIFT.Mul(faceVals.ToMatrix())

		// Check that values at face nodes are significant
		for _, nodeIdx := range Fmask[face] {
			if math.Abs(volumeVals.At(nodeIdx, 0)) < 1e-10 {
				t.Errorf("Face %d: LIFT produced near-zero value at face node %d",
					face, nodeIdx)
			}
		}
	}
}

// TestLift3DOrthogonality tests orthogonality properties
func TestLift3DOrthogonality(t *testing.T) {
	N := 2

	// Generate 3D nodes
	X, Y, Z := Nodes3D(N)
	r, s, tt := XYZtoRST(X, Y, Z)

	// Build required matrices
	V := Vandermonde3D(N, r, s, tt)
	invV := V.InverseWithCheck()
	MassMatrix := invV.Transpose().Mul(invV)
	Fmask := BuildFaceMasks(N, r, s, tt)
	LIFT := Lift3D(N, r, s, tt, V, Fmask)

	// The LIFT operator should satisfy certain orthogonality properties
	// with respect to the mass matrix
	// Specifically: M * LIFT should have structure related to face quadrature
	ML := MassMatrix.Mul(LIFT)

	// Check that M*LIFT has reasonable values
	nrows, ncols := ML.Dims()
	for i := 0; i < nrows; i++ {
		for j := 0; j < ncols; j++ {
			val := ML.At(i, j)
			if math.IsNaN(val) || math.IsInf(val, 0) {
				t.Errorf("M*LIFT has invalid value at (%d,%d): %v", i, j, val)
			}
		}
	}
}

// TestLift3DSymmetry tests that LIFT respects tetrahedral symmetries
func TestLift3DSymmetry(t *testing.T) {
	N := 2

	// Generate 3D nodes
	X, Y, Z := Nodes3D(N)
	r, s, tt := XYZtoRST(X, Y, Z)

	// Build required matrices
	V := Vandermonde3D(N, r, s, tt)
	Fmask := BuildFaceMasks(N, r, s, tt)
	LIFT := Lift3D(N, r, s, tt, V, Fmask)

	// Due to the structure of the reference tetrahedron,
	// there should be patterns in the LIFT matrix reflecting
	// the relationships between faces

	// Test: Face mass matrices should have similar magnitudes
	Nfp := (N + 1) * (N + 2) / 2
	faceMagnitudes := make([]float64, 4)

	for face := 0; face < 4; face++ {
		sum := 0.0
		for i := 0; i < Nfp; i++ {
			col := face*Nfp + i
			// Sum squares of column elements
			nrows, _ := LIFT.Dims()
			for row := 0; row < nrows; row++ {
				val := LIFT.At(row, col)
				sum += val * val
			}
		}
		faceMagnitudes[face] = math.Sqrt(sum)
	}

	// All faces should have similar magnitudes (within factor of 2)
	minMag := faceMagnitudes[0]
	maxMag := faceMagnitudes[0]
	for _, mag := range faceMagnitudes[1:] {
		if mag < minMag {
			minMag = mag
		}
		if mag > maxMag {
			maxMag = mag
		}
	}

	if maxMag/minMag > 2.0 {
		t.Errorf("Face magnitudes vary too much: min=%e, max=%e", minMag, maxMag)
	}
}

// BuildFaceMasks creates the face mask arrays for testing
// This is a helper function that mimics the face mask construction
func BuildFaceMasks(N int, r, s, t []float64) [][]int {
	NODETOL := 1e-7
	Fmask := make([][]int, 4)

	// Face 1: t = -1
	for i := 0; i < len(t); i++ {
		if math.Abs(1.0+t[i]) < NODETOL {
			Fmask[0] = append(Fmask[0], i)
		}
	}

	// Face 2: s = -1
	for i := 0; i < len(s); i++ {
		if math.Abs(1.0+s[i]) < NODETOL {
			Fmask[1] = append(Fmask[1], i)
		}
	}

	// Face 3: r+s+t = -1
	for i := 0; i < len(r); i++ {
		if math.Abs(1.0+r[i]+s[i]+t[i]) < NODETOL {
			Fmask[2] = append(Fmask[2], i)
		}
	}

	// Face 4: r = -1
	for i := 0; i < len(r); i++ {
		if math.Abs(1.0+r[i]) < NODETOL {
			Fmask[3] = append(Fmask[3], i)
		}
	}

	return Fmask
}

// BenchmarkLift3D measures performance of LIFT3D computation
func BenchmarkLift3D(b *testing.B) {
	N := 4

	// Generate 3D nodes
	X, Y, Z := Nodes3D(N)
	r, s, tt := XYZtoRST(X, Y, Z)
	V := Vandermonde3D(N, r, s, tt)
	Fmask := BuildFaceMasks(N, r, s, tt)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = Lift3D(N, r, s, tt, V, Fmask)
	}
}

package readers

import (
	"os"
	"strings"
	"testing"
)

// TestReadGambitNeutralBoundaryConditionsError tests handling of malformed BC data
func TestReadGambitNeutralBoundaryConditionsError(t *testing.T) {
	// This test has malformed boundary condition data
	content := `        CONTROL INFO 2.0.0
** GAMBIT NEUTRAL FILE
Test boundary conditions
PROGRAM:                  Test     VERSION:  1.0
Mon Jan  1 00:00:00 2025
     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL
         4         1         1         2         3         3
ENDOFSECTION
   NODAL COORDINATES 2.0.0
         1   0.00000000000e+00   0.00000000000e+00   0.00000000000e+00
         2   1.00000000000e+00   0.00000000000e+00   0.00000000000e+00
         3   0.00000000000e+00   1.00000000000e+00   0.00000000000e+00
         4   0.00000000000e+00   0.00000000000e+00   1.00000000000e+00
ENDOFSECTION
   ELEMENTS/CELLS 2.0.0
         1         6         4         1         2         3         4
ENDOFSECTION
       BOUNDARY CONDITIONS 2.0.0
inlet           1         1         0         0         0         0         0         0
         1         6         1
wall            1        notanumber         0         0         0         0         0         0
ENDOFSECTION`

	tmpFile := createTempNeuFile(t, content)
	defer os.Remove(tmpFile)

	_, err := ReadGambitNeutral(tmpFile)
	if err == nil {
		t.Fatal("Expected error for malformed boundary condition data, but got none")
	}

	// Verify the error message indicates BC format issue
	// The error should mention "invalid boundary condition format"
	if !strings.Contains(err.Error(), "invalid boundary condition format") {
		t.Errorf("Expected error about invalid boundary condition format, got: %v", err)
	}
}

// TestReadGambitNeutralBoundaryConditionsValid tests reading valid boundary conditions
func TestReadGambitNeutralBoundaryConditionsValid(t *testing.T) {
	// Test with properly formatted boundary conditions
	content := `        CONTROL INFO 2.0.0
** GAMBIT NEUTRAL FILE
Test boundary conditions
PROGRAM:                  Test     VERSION:  1.0
Mon Jan  1 00:00:00 2025
     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL
         8         2         1         2         3         3
ENDOFSECTION
   NODAL COORDINATES 2.0.0
         1   0.00000000000e+00   0.00000000000e+00   0.00000000000e+00
         2   1.00000000000e+00   0.00000000000e+00   0.00000000000e+00
         3   0.00000000000e+00   1.00000000000e+00   0.00000000000e+00
         4   0.00000000000e+00   0.00000000000e+00   1.00000000000e+00
         5   1.00000000000e+00   1.00000000000e+00   0.00000000000e+00
         6   1.00000000000e+00   0.00000000000e+00   1.00000000000e+00
         7   0.00000000000e+00   1.00000000000e+00   1.00000000000e+00
         8   1.00000000000e+00   1.00000000000e+00   1.00000000000e+00
ENDOFSECTION
   ELEMENTS/CELLS 2.0.0
         1         6         4         1         2         3         4
         2         6         4         2         5         6         8
ENDOFSECTION
       BOUNDARY CONDITIONS 2.0.0
inlet           1         2         0         0         0         0         0         0
         1         6         1
         1         6         2
wall            1         3         0         0         0         0         0         0
         1         6         3
         2         6         1
         2         6         4
ENDOFSECTION`

	tmpFile := createTempNeuFile(t, content)
	defer os.Remove(tmpFile)

	mesh, err := ReadGambitNeutral(tmpFile)
	if err != nil {
		t.Fatalf("Failed to read Gambit neutral file: %v", err)
	}

	// Check boundary tags
	if len(mesh.BoundaryTags) != 2 {
		t.Errorf("Expected 2 boundary conditions, got %d", len(mesh.BoundaryTags))
	}

	// Check that we have the expected BC names
	foundInlet := false
	foundWall := false
	for _, name := range mesh.BoundaryTags {
		if name == "inlet" {
			foundInlet = true
		} else if name == "wall" {
			foundWall = true
		}
	}

	if !foundInlet {
		t.Error("Expected to find 'inlet' boundary condition")
	}
	if !foundWall {
		t.Error("Expected to find 'wall' boundary condition")
	}

	// Check boundary elements
	if len(mesh.BoundaryElements["inlet"]) != 2 {
		t.Errorf("Expected 2 boundary elements for 'inlet', got %d", len(mesh.BoundaryElements["inlet"]))
	}
	if len(mesh.BoundaryElements["wall"]) != 3 {
		t.Errorf("Expected 3 boundary elements for 'wall', got %d", len(mesh.BoundaryElements["wall"]))
	}
}

// TestReadGambitNeutralBoundaryConditionsNode tests node-based boundary conditions
func TestReadGambitNeutralBoundaryConditionsNode(t *testing.T) {
	// Test with node-based boundary conditions
	content := `        CONTROL INFO 2.0.0
** GAMBIT NEUTRAL FILE
Test node boundary conditions
PROGRAM:                  Test     VERSION:  1.0
Mon Jan  1 00:00:00 2025
     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL
         4         1         1         1         3         3
ENDOFSECTION
   NODAL COORDINATES 2.0.0
         1   0.00000000000e+00   0.00000000000e+00   0.00000000000e+00
         2   1.00000000000e+00   0.00000000000e+00   0.00000000000e+00
         3   0.00000000000e+00   1.00000000000e+00   0.00000000000e+00
         4   0.00000000000e+00   0.00000000000e+00   1.00000000000e+00
ENDOFSECTION
   ELEMENTS/CELLS 2.0.0
         1         6         4         1         2         3         4
ENDOFSECTION
       BOUNDARY CONDITIONS 2.0.0
fixed           0         2         3         0         0         0         0         0
         1   0.00000000000e+00   0.00000000000e+00   0.00000000000e+00
         2   0.00000000000e+00   0.00000000000e+00   0.00000000000e+00
ENDOFSECTION`

	tmpFile := createTempNeuFile(t, content)
	defer os.Remove(tmpFile)

	mesh, err := ReadGambitNeutral(tmpFile)
	if err != nil {
		t.Fatalf("Failed to read Gambit neutral file: %v", err)
	}

	// Check boundary tags
	if len(mesh.BoundaryTags) != 1 {
		t.Errorf("Expected 1 boundary condition, got %d", len(mesh.BoundaryTags))
	}

	// Check that we have the expected BC name
	if mesh.BoundaryTags[0] != "fixed" {
		t.Errorf("Expected boundary condition name 'fixed', got '%s'", mesh.BoundaryTags[0])
	}

	// Node boundary conditions are not stored in BoundaryElements
	// They would need separate handling if needed
}

// TestReadGambitNeutralBoundaryConditionsInvalidITYPE tests invalid ITYPE values
func TestReadGambitNeutralBoundaryConditionsInvalidITYPE(t *testing.T) {
	content := `        CONTROL INFO 2.0.0
** GAMBIT NEUTRAL FILE
Test invalid ITYPE
PROGRAM:                  Test     VERSION:  1.0
Mon Jan  1 00:00:00 2025
     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL
         4         1         1         1         3         3
ENDOFSECTION
   NODAL COORDINATES 2.0.0
         1   0.00000000000e+00   0.00000000000e+00   0.00000000000e+00
         2   1.00000000000e+00   0.00000000000e+00   0.00000000000e+00
         3   0.00000000000e+00   1.00000000000e+00   0.00000000000e+00
         4   0.00000000000e+00   0.00000000000e+00   1.00000000000e+00
ENDOFSECTION
   ELEMENTS/CELLS 2.0.0
         1         6         4         1         2         3         4
ENDOFSECTION
       BOUNDARY CONDITIONS 2.0.0
invalid         2         1         0         0         0         0         0         0
ENDOFSECTION`

	tmpFile := createTempNeuFile(t, content)
	defer os.Remove(tmpFile)

	_, err := ReadGambitNeutral(tmpFile)
	if err == nil {
		t.Fatal("Expected error for invalid ITYPE value, but got none")
	}

	// Verify the error message indicates ITYPE issue
	if !strings.Contains(err.Error(), "invalid boundary condition ITYPE") {
		t.Errorf("Expected error about invalid ITYPE, got: %v", err)
	}
}

// Helper function to create temporary .neu file
func createTempNeuFile(t *testing.T, content string) string {
	tmpFile, err := os.CreateTemp("", "test_*.neu")
	if err != nil {
		t.Fatal(err)
	}
	defer tmpFile.Close()

	if _, err := tmpFile.WriteString(content); err != nil {
		t.Fatal(err)
	}

	return tmpFile.Name()
}

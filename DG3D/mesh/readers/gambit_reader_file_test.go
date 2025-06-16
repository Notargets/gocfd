package readers

import (
	"fmt"
	"os"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

func TestReadGambitNeutral_WithPartitions(t *testing.T) {
	// Test reading the cube-partitioned.neu file
	meshPath := "../cube-partitioned.neu" // Adjust path as needed

	mesh, err := ReadGambitNeutral(meshPath)
	require.NoError(t, err, "Failed to read partitioned mesh")
	require.NotNil(t, mesh)

	// Verify EToP was populated
	require.NotNil(t, mesh.EToP, "EToP should be populated for partitioned mesh")
	require.Equal(t, mesh.NumElements, len(mesh.EToP), "EToP should have entry for each element")

	// Count elements per partition
	partitionCounts := make(map[int]int)
	maxPartID := -1
	for _, p := range mesh.EToP {
		if p >= 0 {
			partitionCounts[p]++
			if p > maxPartID {
				maxPartID = p
			}
		}
	}

	// Should have 4 partitions (p1, p2, p3, p4) based on the file
	assert.Equal(t, 4, len(partitionCounts), "Should have 4 partitions")
	assert.Equal(t, 4, maxPartID, "Max partition ID should be 4")

	// Each partition should have roughly equal number of elements
	// File shows: p1=142, p2=141, p3=141, p4=141 elements
	for p, count := range partitionCounts {
		t.Logf("Partition %d: %d elements", p, count)
		assert.Greater(t, count, 100, "Each partition should have > 100 elements")
		assert.Less(t, count, 150, "Each partition should have < 150 elements")
	}

	// Verify total element count
	totalPartitioned := 0
	for _, count := range partitionCounts {
		totalPartitioned += count
	}
	assert.Equal(t, mesh.NumElements, totalPartitioned, "All elements should be partitioned")

	// Verify element groups contain partition info
	foundPartitionGroups := 0
	for _, group := range mesh.ElementGroups {
		if len(group.Name) >= 2 && group.Name[0] == 'p' {
			// This is a partition group
			foundPartitionGroups++

			// Extract expected partition ID
			partID := extractPartitionIDFromName(group.Name)
			if partID >= 0 {
				// Verify all elements in this group have correct partition assignment
				for _, elemIdx := range group.Elements {
					assert.Equal(t, partID, mesh.EToP[elemIdx],
						"Element %d in group %s should be in partition %d",
						elemIdx, group.Name, partID)
				}
			}
		}
	}
	assert.Equal(t, 4, foundPartitionGroups, "Should find 4 partition groups (p1-p4)")
}

func TestReadGambitNeutral_NonPartitionedMesh(t *testing.T) {
	// Create a simple non-partitioned mesh for comparison
	content := `        CONTROL INFO 2.0.0
** GAMBIT NEUTRAL FILE
Test mesh without partitions
PROGRAM:                  Test     VERSION:  1.0
Mon Jan  1 00:00:00 2025
     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL
         4         1         1         0         3         3
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
       ELEMENT GROUP 2.0.0
GROUP:          1 ELEMENTS:          1 MATERIAL:          0 NFLAGS:          1
                Material group 1
       0
       1
ENDOFSECTION`

	tmpFile := createTempNeuFile(t, content)
	defer os.Remove(tmpFile)

	mesh, err := ReadGambitNeutral(tmpFile)
	require.NoError(t, err)

	// For non-partitioned mesh (single group), EToP should be nil
	assert.Nil(t, mesh.EToP, "Non-partitioned mesh should have nil EToP")
}

func TestReadGambitNeutral_PartitionedMeshStructure(t *testing.T) {
	// Test with a small partitioned mesh to verify structure
	content := `        CONTROL INFO 2.0.0
** GAMBIT NEUTRAL FILE
Test partitioned mesh
PROGRAM:                  Test     VERSION:  1.0
Mon Jan  1 00:00:00 2025
     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL
         8         4         2         0         3         3
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
         2         6         4         2         5         3         6
         3         6         4         3         4         7         8
         4         6         4         5         6         7         8
ENDOFSECTION
       ELEMENT GROUP 2.0.0
GROUP:          1 ELEMENTS:          2 MATERIAL:          0 NFLAGS:          1
                              p1
       0
       1         2
ENDOFSECTION
       ELEMENT GROUP 2.0.0
GROUP:          2 ELEMENTS:          2 MATERIAL:          0 NFLAGS:          1
                              p2
       0
       3         4
ENDOFSECTION`

	tmpFile := createTempNeuFile(t, content)
	defer os.Remove(tmpFile)

	mesh, err := ReadGambitNeutral(tmpFile)
	require.NoError(t, err)

	// Verify partitions were detected
	require.NotNil(t, mesh.EToP, "EToP should be populated for 2-group mesh")
	require.Equal(t, 4, len(mesh.EToP), "Should have 4 elements")

	// Check partition assignments
	assert.Equal(t, 1, mesh.EToP[0], "Element 0 should be in partition 1")
	assert.Equal(t, 1, mesh.EToP[1], "Element 1 should be in partition 1")
	assert.Equal(t, 2, mesh.EToP[2], "Element 2 should be in partition 2")
	assert.Equal(t, 2, mesh.EToP[3], "Element 3 should be in partition 2")
}

// Helper function to extract partition ID from group name
func extractPartitionIDFromName(name string) int {
	if len(name) >= 2 && name[0] == 'p' {
		// Try to parse the number after 'p'
		var partID int
		_, err := fmt.Sscanf(name[1:], "%d", &partID)
		if err == nil {
			return partID
		}
	}
	return -1
}

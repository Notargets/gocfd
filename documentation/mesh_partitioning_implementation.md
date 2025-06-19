# Mesh Partitioning Design and Implementation Guide

## Table of Contents

1. [Overview](#overview)
2. [Mesh Splitting Architecture](#mesh-splitting-architecture)
3. [Face Buffer Mechanism](#face-buffer-mechanism)
4. [Data Structures](#data-structures)
5. [Implementation Details](#implementation-details)
6. [Test Infrastructure](#test-infrastructure)
7. [Unit Testing Strategy](#unit-testing-strategy)

## Overview

The mesh partitioning system splits a global mesh into multiple partitions for
parallel computation while maintaining:

- Element connectivity across partition boundaries
- Boundary condition mappings
- Geometric fidelity
- Face communication patterns

### Key Design Principles

- **Separation of Concerns**: Mesh splitting handles geometry; face buffers
  handle inter-partition communication
- **Minimal Redundancy**: Reuse existing DG3D initialization for each partition
- **Testability**: Progressive validation from simple to complex cases

## Mesh Splitting Architecture

### Workflow

```
Global Mesh → Partition Assignment (EToP) → Split Mesh → Local DG3D Initialization
                                                ↓
                                         Face Buffer Setup ← EToE + EToP
```

### Key Components

#### 1. MeshSplitter

Responsible for dividing the global mesh into partition-local meshes.

```go
type MeshSplitter struct {
    // Inputs
    GlobalMesh *mesh.Mesh
    VX, VY, VZ []float64 // Global vertex coordinates
    EToP       []int     // Element to partition mapping

    // Outputs
    PartitionMeshes []mesh.Mesh
    PEToE          map[int][]int // Partition element to global element
    EToP_0based    []int             // Normalized to 0-based indexing

    // Internal mappings
    vertexMaps map[int]map[int]int // partID -> globalVertexID -> localVertexID
}
```

#### 2. Partition Mesh Generation Process

1. **Normalize EToP**: Convert to 0-based indexing if needed
2. **Build PEToE**: Create inverse mapping for validation
3. **Extract Elements**: Group elements by partition
4. **Extract Vertices**: Collect unique vertices per partition
5. **Create Local Mesh**: Build mesh.Mesh with local indexing

## Face Buffer Mechanism

The face buffer system efficiently manages face data exchange between partitions
without explicit index matching.

### Build Phase (FaceBufferBuilder)

```go
// Input: Element connectivity + partition assignment
EToE [][]int // Element-to-element connectivity
EToP []int    // Element-to-partition mapping

// Automatic face classification:
for each element E:
    for each face F:
        if EToE[E][F] == E:
            → Boundary face else if EToP[E] == EToP[EToE[E][F]]:
            → Local neighbor (same partition) else:
            → Remote neighbor (different partition)
```

### Runtime Phase (FaceBufferRuntime)

The runtime structure provides efficient traversal of face points:

```go
// Natural traversal order: element → face → point
for elem := 0; elem < K; elem++ {
    for face := 0; face < Nface; face++ {
        for point := 0; point < Nfp; point++ {
            mValue := GetMValue(equation)
            pValue := GetPValue(equation) // Automatically handles all face types
// Apply Riemann solver or boundary condition
} } }
```

### Communication Pattern

```
Partition A                     Partition B
SendBuffer[indices] ----MPI----> ReceiveBuffer
    ↑                                ↓
Assembly from                   Used by GetPValue()
local M array                   for remote faces
```

## Data Structures

### Core Mesh Structures

```go
// Global mesh representation
type Mesh struct {
    // Geometry
    Vertices     [][]float64
    Elements     [][]int
    ElementTypes []ElementType

    // Connectivity
    EToE [][]int // Element to element
    EToF [][]int // Element to face
    EToP []int   // Element to partition

    // Boundary conditions
    BoundaryElements map[string][][2]int // bcTag -> [][elemID, faceID]
}

// Partition-specific BC mapping
type PartitionBCMap struct {
    BCTag    string
    Elements []int // Local element IDs
    Faces    []int // Local face IDs
}
```

### Face Buffer Structures

```go
// Build-time connection information
type BuildTimeConnection struct {
    MPos        uint32   // M-side position
    Type        FaceType // Boundary/Local/Remote
    BCType      uint32   // BC type if boundary
    PartitionID uint32 // Remote partition if remote
    PPos        uint32 // P-side position
}

// Runtime arrays for efficient execution
type FaceBufferRuntime struct {
    // Flat arrays indexed by M position
    FaceTypes    []FaceType
    BCTypes      []uint32
    PartitionIDs []uint32
    LocalPIndices []uint32

    // Remote communication buffers
    RemoteBuffers map[uint32]*RemoteBufferData
}
```

## Implementation Details

### SplitMesh Method

```go
func (ms *MeshSplitter) SplitMesh() ([]mesh.Mesh, map[int][]int, error) {
    // Step 1: Normalize EToP to 0-based
    ms.normalizeEToP()

    // Step 2: Build PEToE mapping
    ms.buildPEToE()

    // Step 3: Create partition meshes
    partitionMeshes := make([]mesh.Mesh, ms.numPartitions)

    for partID := 0; partID < ms.numPartitions; partID++ {
        // Extract elements for this partition
        globalElemIDs := ms.PEToE[partID]

        // Extract required vertices
        localVX, localVY, localVZ, vertexMap := ms.extractPartitionVertices(
        partID, globalElemIDs)

        // Build local element connectivity
        localElements := ms.buildLocalElements(globalElemIDs, vertexMap)

        // Create partition mesh
        partitionMeshes[partID] = ms.createPartitionMesh(
        localVX, localVY, localVZ, localElements)

        // Transform BC mappings
        ms.transformBCMappings(partID, globalElemIDs)
    }

    return partitionMeshes, ms.PEToE, nil
}
```

### BC Transformation

```go
func (ms *MeshSplitter) transformBCMappings(partID int, globalElemIDs []int) {
    // Create global to local element mapping
    globalToLocal := make(map[int]int)
    for localIdx, globalIdx := range globalElemIDs {
        globalToLocal[globalIdx] = localIdx
    }

    // Transform each BC
    for bcTag, globalBCs := range ms.GlobalMesh.BoundaryElements {
        for _, bc := range globalBCs {
            globalElemID, faceID := bc[0], bc[1]

            if localElemID, exists := globalToLocal[globalElemID]; exists {
                // This BC belongs to this partition
                ms.addPartitionBC(partID, bcTag, localElemID, faceID)
            }
        }
    }
}
```

## Test Infrastructure

### Enhanced Test Mesh Helpers

Extend `utils/mesh_test_helpers.go` with partitioning and BC support:

```go
// TestMeshConfig provides configuration for test mesh generation
type TestMeshConfig struct {
    MeshType      string // "SingleTet", "TwoTets", "CubeMesh", etc.
    PartitionType string // "None", "Bisect", "ByElement", "Custom"
    Partitions    []int               // Custom partition assignment
    BCMappings    map[string][][2]int // BC tag -> [][elemID, faceID]
}

// StandardTestCase provides reusable test configurations
type StandardTestCase struct {
    Name          string
    Config        TestMeshConfig
    ExpectedProps TestProperties
}

// TestProperties for validation
type TestProperties struct {
    NumPartitions      int
    ElementsPerPart    []int
    SharedFaces        int
    BoundaryFaces      int
    InterPartitionComm map[[2]int]int // [partA,partB] -> num shared faces
}

// GetStandardTestCases returns predefined test cases
func GetStandardTestCases() []StandardTestCase {
	return []StandardTestCase{
		{
			Name: "SingleTetSinglePartition",
			Config: TestMeshConfig{
				MeshType:      "SingleTet",
				PartitionType: "None",
				BCMappings: map[string][][2]int{
					"Wall": {{0, 0}, {0, 1}, {0, 2}, {0, 3}},
				},
			},
			ExpectedProps: TestProperties{
				NumPartitions:   1,
				ElementsPerPart: []int{1},
				BoundaryFaces:   4,
			},
		},
		{
			Name: "TwoTetsBisected",
			Config: TestMeshConfig{
				MeshType:      "TwoTets",
				PartitionType: "Bisect",
				BCMappings: map[string][][2]int{
					"Inflow":  {{0, 0}},
					"Outflow": {{1, 3}},
				},
			},
			ExpectedProps: TestProperties{
				NumPartitions:      2,
				ElementsPerPart:    []int{1, 1},
				SharedFaces:        1,
				BoundaryFaces:      6,
				InterPartitionComm: map[[2]int]int{{0, 1}: 1},
			},
		},
		{
			Name: "CubeMeshQuadPartition",
			Config: TestMeshConfig{
				MeshType:      "CubeMesh",
				PartitionType: "Custom",
				Partitions:    []int{0, 0, 1, 1, 2, 3}, // 6 tets → 4 partitions
			},
			ExpectedProps: TestProperties{
				NumPartitions:   4,
				ElementsPerPart: []int{2, 2, 1, 1},
			},
		},
	}
}
```

### Test Mesh Factory

```go
type TestMeshFactory struct {
    baseHelpers *TestMeshes
}

func (tmf *TestMeshFactory) CreateTestMesh(config TestMeshConfig) (*mesh.Mesh, error) {
    // Get base mesh
    var completeMesh CompleteMesh
    switch config.MeshType {
        case "SingleTet":
            completeMeshtype TestMeshFactory struct {
	            baseHelpers *TestMeshes
            }

            func (tmf *TestMeshFactory) CreateTestMesh(config TestMeshConfig) (*mesh.Mesh, error) {
	            // Get base mesh
	            var completeMesh CompleteMesh
	            switch config.MeshType {
	            case "SingleTet":
		            completeMesh = tmf.createSingleTetMesh()
	            case "TwoTets":
		            completeMesh = tmf.baseHelpers.TwoTetMesh
	            case "CubeMesh":
		            completeMesh = tmf.baseHelpers.CubeMesh
		            // ... other cases
	            }

	            // Convert to mesh.Mesh
	            m := ConvertToMesh(completeMesh)

	            // Apply partitioning
	            m.EToP = tmf.applyPartitioning(m, config.PartitionType, config.Partitions)

	            // Apply boundary conditions
	            m.BoundaryElements = config.BCMappings

	            return m, nil
            } = tmf.createSingleTetMesh()
        case "TwoTets":
            completeMesh = tmf.baseHelpers.TwoTetMesh
        case "CubeMesh":
            completeMesh = tmf.baseHelpers.CubeMesh
        // ... other cases
        }

    // Convert to mesh.Mesh
    m := ConvertToMesh(completeMesh)

    // Apply partitioning
    m.EToP = tmf.applyPartitioning(m, config.PartitionType, config.Partitions)

    // Apply boundary conditions
    m.BoundaryElements = config.BCMappings

    return m, nil
}
```

## Unit Testing Strategy

Following the Unit Testing Principles, tests progress from fundamental to
complex:

### Level 1: Foundation Tests

```go
// Test EToP normalization
func TestEtoPNormalization(t *testing.T) {
    testCases := []struct {
        name     string
        input    []int
        expected []int
    }{
        {"AlreadyZeroBased", []int{0, 1, 2}, []int{0, 1, 2}},
        {"OneBased", []int{1, 2, 3}, []int{0, 1, 2}},
        {"Arbitrary", []int{5, 7, 5, 9}, []int{0, 2, 0, 4}},
    }
    // ... test implementation
}
```

### Level 2: Single Partition Validation

```go
// Test identity operation: single partition should preserve mesh
func TestSinglePartitionIdentity(t *testing.T) {
    tc := GetStandardTestCases()[0] // SingleTetSinglePartition

    // Split mesh
    splitter := NewMeshSplitter(globalMesh, vx, vy, vz, eToP)
    partMeshes, peToE, _ := splitter.SplitMesh()

    // Verify identity preservation
    assert.Equal(t, 1, len(partMeshes))
    assert.Equal(t, globalMesh.NumElements, partMeshes[0].NumElements)
    assert.Equal(t, globalMesh.NumVertices, partMeshes[0].NumVertices)
}
```

### Level 3: Partition Connectivity

```go
// Test face buffer connectivity between partitions
func TestFaceBufferConnectivity(t *testing.T) {
    tc := GetStandardTestCases()[1] // TwoTetsBisected

    // Create face buffers with test data
    for partID, partMesh := range partitionMeshes {
        fb := NewFaceBufferBuilder(...)

        // Load with global face IDs for testing
        for elem := 0; elem < partMesh.NumElements; elem++ {
            globalElem := peToE[partID][elem]
            for face := 0; face < 4; face++ {
                for point := 0; point < nfp; point++ {
                    globalFaceID := computeGlobalFaceID(globalElem, face)
                    fb.SetMValue(elem, face, point, 0, float32(globalFaceID))
                }
            }
        }
    }

    // Exchange buffers and verify
    // ... communication simulation and validation
}
```

### Level 4: Mathematical Properties

```go
// Test volume conservation across partitions
func TestVolumeConservation(t *testing.T) {
    for _, tc := range GetStandardTestCases() {
        t.Run(tc.Name, func (t *testing.T) {
            // Compute global mesh volume
            globalVolume := computeMeshVolume(globalMesh)

            // Compute sum of partition volumes
            partitionVolume := 0.0
            for _, partMesh := range partitionMeshes {
                partitionVolume += computeMeshVolume(partMesh)
            }

            // Verify conservation
            assert.InDelta(t, globalVolume, partitionVolume, 1e-12)
        })
	}
}
```

### Level 5: BC Mapping Correctness

```go
// Test BC transformation from global to local coordinates
func TestBCMapping(t *testing.T) {
    // For each BC in global mesh
    for bcTag, globalBCs := range globalMesh.BoundaryElements {
        for _, bc := range globalBCs {
            globalElem, globalFace := bc[0], bc[1]

            // Find partition containing this element
            partID := eToP[globalElem]
            localElem := findLocalElement(peToE[partID], globalElem)

            // Verify BC exists in partition
            found := false
            for _, localBC := range partitionBCMaps[partID][bcTag] {
                if localBC[0] == localElem && localBC[1] == globalFace {
                    found = true
                    break
                }
            }
            assert.True(t, found, "BC not found in partition")
        }
    }
}
```

### Progressive Test Structure

1. **Unit Tests**: Test individual functions (normalization, mapping)
2. **Integration Tests**: Test mesh splitting with simple cases
3. **System Tests**: Test face buffer communication
4. **Validation Tests**: Test mathematical properties
5. **Stress Tests**: Test with complex meshes and many partitions

### Test Helpers

```go
// Validation utilities
func validatePartitionCompleteness(global *mesh.Mesh, parts []mesh.Mesh, peToE map[int][]int) error
func validateFaceConnectivity(parts []mesh.Mesh, faceBuffers []*FaceBufferRuntime) error
func validateBCCompleteness(globalBC, partBC map[string][][2]int, peToE map[int][]int) error

// Computation utilities  
func computeGlobalFaceID(elemID, faceID int) int
func computeMeshVolume(m *mesh.Mesh) float64
func computeElementVolume(vertices [][]float64, connectivity []int) float64
```

## Summary

The mesh partitioning system achieves efficient parallel computation through:

1. **Clean separation**: Mesh splitting handles geometry; face buffers handle
   communication
2. **Automatic classification**: Face buffer builder automatically identifies
   partition boundaries
3. **Preserved correctness**: BC mappings and geometric properties maintained
4. **Testable design**: Progressive validation from simple to complex cases
5. **Reusable infrastructure**: Standard test meshes minimize duplication

The design prioritizes correctness and testability while maintaining efficiency
for parallel execution.
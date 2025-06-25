# Face Buffer Design and Implementation

## Overview

The face buffer system provides a solution for managing face data in
partitioned meshes, enabling efficient computation of numerical fluxes through
Riemann solvers. The design achieves high performance by organizing data access
patterns to match the natural flow of computation, eliminating the need for
complex indexing during runtime.

## Implementation Strategy

The face buffer system is built in two phases:

### Phase 1: Face Connectivity Building
- Identifies face types (boundary/interior/remote)
- Builds connectivity indices without specific boundary condition information
- Uses placeholder values for boundary faces

### Phase 2: Boundary Condition Overlay
- Replaces boundary placeholders with specific BC codes
- Performed as a separate pass after connectivity is established
- Allows BC updates without rebuilding connectivity

## Face Buffer Structure and Natural Ordering

### The M Array

At the heart of the face buffer system is the M array, which stores all face
point data for a partition in a specific order called "natural traversal order."
This ordering follows the nested loop structure:

- Outermost loop: Elements (0 to K-1)
- Middle loop: Faces per element (0 to Nfaces-1)
- Innermost loop: Points per face (depends on polynomial order)

This creates a linear sequence where face points are accessed in a predictable
pattern. For example, in a mesh with 100 elements, Nfaces faces per element, and Nfp
points per face, the M array contains 100×Nfaces×Nfp face points arranged sequentially.

### Face Classification

During the build phase, each face (not individual face points) is classified as one of
three types:

1. **Boundary faces**: Faces on the domain boundary requiring boundary conditions
2. **Interior faces**: Faces connecting to other elements within the same partition
3. **Remote faces**: Faces connecting to elements in different partitions

This classification is encoded in a compact face index array of dimension (Nfaces × K_part).

## Index Structures for Face Connectivity

### The Face Index Array

A single integer array of dimension (Nfaces × K_part) encodes all face connectivity information:

```
FaceIndex[face + elem*Nfaces] = encoding
```

Where the encoding uses:
- **Positive values**: Index to the start of P values in the local M buffer (interior faces)
- **Value -999**: Boundary face placeholder (replaced with specific BC codes in Phase 2)
- **Negative values -1 to -998**: Specific boundary condition codes (assigned in Phase 2)
- **Special value -9999**: Remote face indicator

This compact representation eliminates branching during traversal while providing all necessary connectivity information.

### Building the Face Index (Phase 1)

During the initial build phase, the face buffer builder:
1. Analyzes mesh connectivity (VmapM/VmapP, EToE, EToP)
2. For each face, determines type using:
   - If `EToE[elem][face] == elem`: Boundary face
   - If `EToP[EToE[elem][face]] != EToP[elem]`: Remote face
   - Otherwise: Interior face
3. Assigns initial values:
   - Interior faces: Positive index pointing to neighbor's face start in M buffer
   - Boundary faces: -999 (placeholder)
   - Remote faces: -9999

Note: EToP array must contain partition assignments for all global elements, not just local ones.

### Boundary Condition Overlay (Phase 2)

After connectivity is established, a separate pass:
1. Scans FaceIndex for -999 values
2. Uses mesh boundary condition data to determine specific BC type
3. Replaces -999 with actual BC code (-1 through -998)

This two-phase approach separates connectivity concerns from boundary condition specifics.

### LocalPIndices for Interior Faces

For faces with positive index values:
- The value points to the START of the Nfp P values in the same partition's M buffer
- Example: If FaceIndex[f] = 1200, then P values for point i are at M[1200 + i]
- Built during preprocessing by finding where the neighbor element's matching face appears in the M buffer
- The M buffer serves dual purpose: source of M values during traversal AND source of P values for interior connections

### RemoteSendIndices for Cross-Partition Communication

Each partition maintains a map structure:
```
RemoteSendIndices[targetPartition] = [list of M buffer positions]
```

This structure:
- Identifies which face values from THIS partition's M buffer are needed by OTHER partitions
- Built during preprocessing by analyzing cross-partition face connections
- Used to pack send buffers before each time step
- Ensures values are packed in the exact order the receiving partition will consume them

## Runtime Traversal for Riemann Solvers

### The Fundamental Operation

When computing numerical fluxes, the Riemann solver requires two values at each
face point:

- The "minus" value (M): from the current element
- The "plus" value (P): from the neighboring element

The face buffer system retrieves these values efficiently during a single linear
traversal with minimal branching.

### Sequential Access with Face-Level Indexing

The runtime traversal follows this pattern:

```
for elem = 0 to K_part-1:
    for face = 0 to Nfaces-1:
        face_code = FaceIndex[face + elem*Nfaces]
        
        for point = 0 to Nfp-1:
            m_idx = elem*Nfaces*Nfp + face*Nfp + point
            M = M_buffer[m_idx]
            
            if face_code > 0:  // Interior face
                P = M_buffer[face_code + point]
            elif face_code == -9999:  // Remote face
                P = receive_buffer[partition][remote_counter++]
            else:  // Boundary face (BC code)
                P = applyBC(M, -face_code)
            
            flux[m_idx] = RiemannSolver(M, P)
```

This design ensures:
- Sequential memory access through the M buffer
- Minimal branching (one check per face, not per point)
- Direct indexing for interior P values
- Sequential consumption of remote receive buffers

## Communication Pattern for Remote Faces

### Pre-traversal Communication

Before the Riemann solver traversal begins:

1. **Pack Send Buffers**: Each partition uses RemoteSendIndices to identify which M values to send:
   ```
   for each targetPartition:
       for each idx in RemoteSendIndices[targetPartition]:
           sendBuffer[targetPartition].append(M_buffer[idx])
   ```

2. **Exchange Buffers**: Send buffers are exchanged between partitions (typically using MPI)

3. **Receive Buffers Ready**: Each partition now has receive buffers containing P values from other partitions, ordered for sequential consumption

### Natural Order Preservation

The key insight is that RemoteSendIndices are ordered to match the receiving partition's traversal order:
- No reordering needed during packing or unpacking
- Receive buffers are consumed sequentially with a simple counter per partition
- The sender pre-sorts the data in the order the receiver needs it

## Boundary Processing in Partitioned Meshes

### Efficient BC Encoding

Boundary conditions are encoded as negative values in the face index:
- BC type 1: -1
- BC type 2: -2
- ...
- BC type 998: -998
- Boundary placeholder: -999 (Phase 1 only)
- Remote face: -9999

This allows:
- Single integer comparison to determine face type
- Direct BC type extraction: bc_type = -face_code
- No separate arrays needed for face classification

### Boundary Condition Application

When a boundary face is encountered (negative face_code ∈ [-998, -1]):
1. Extract BC type from the negative code
2. Apply appropriate formula based on BC type
3. Common examples:
   - Wall (BC 1): P = -M (reflection)
   - Outflow (BC 2): P = M (zero gradient)
   - Inflow (BC 3): P = prescribed value

The placeholder value -999 should never appear during runtime traversal, as it's replaced during the BC overlay phase.

## Performance Benefits

This design achieves several performance advantages:

1. **Compact Storage**: Single integer array encodes all face connectivity
2. **Sequential Access**: Both M array and receive buffers accessed linearly
3. **Minimal Branching**: One branch per face (not per face point)
4. **Cache Efficiency**: Predictable access patterns and compact data structures
5. **Direct Indexing**: No searching or indirect lookups during flux computation

## Example: Processing Element 15

Consider processing element 15 with Nfaces=4 faces in a partition with Nfp=6:

```
Face 0: FaceIndex[15*4 + 0] = 720    // Interior to element 30, face 0
Face 1: FaceIndex[15*4 + 1] = -1     // Wall boundary
Face 2: FaceIndex[15*4 + 2] = -9999  // Remote connection to partition 2
Face 3: FaceIndex[15*4 + 3] = 1680   // Interior to element 70, face 2
```

During traversal:
- Face 0: P values from M[720] through M[725]
- Face 1: P values computed from wall BC formula
- Face 2: P values from receive_buffer[2][counter_2] through [counter_2+5]
- Face 3: P values from M[1680] through M[1685]

## Implementation Details

### Data Structure Evolution

The face buffer builder creates:
```go
type FaceBuffer struct {
    // Dimensions
    Nfp      uint32  // Face points per face
    Nfaces   uint32  // Faces per element
    K        uint32  // Number of elements in this partition
    
    // Face connectivity (Phase 1)
    FaceIndex []int32  // Size: Nfaces × K, face-level connectivity
    
    // Remote communication
    RemoteSendIndices map[uint32][]uint32  // [partitionID] -> M buffer indices
}
```

Key features:
- FaceIndex has dimension (Nfaces × K), not (Nfp × Nfaces × K)
- No separate LocalPIndices array - embedded as positive values in FaceIndex
- Each index points to the start of Nfp contiguous values

### Helper Functions

```go
// Phase 1: Build connectivity
func BuildFaceBuffer(el *Element3D) (*FaceBuffer, error)

// Phase 2: Apply BC overlay
func ApplyBoundaryConditions(fb *FaceBuffer, bcData BCData) error
```

### Index Calculation Example

For element e, face f, point p:
```
// M value location (always direct)
m_idx = e*Nfaces*Nfp + f*Nfp + p

// P value location (depends on face type)
face_code = FaceIndex[f + e*Nfaces]
if face_code > 0:
    p_idx = face_code + p  // Interior: offset from face start
else if face_code == -9999:
    p_idx = next from receive buffer  // Remote
else:
    p = BC formula using code  // Boundary
```

## Summary

The face buffer system transforms complex face data management into efficient sequential operations through:

1. **Two-phase construction**: Connectivity building followed by BC overlay
2. **Dual-purpose M buffer**: Serves as both M source and P source for interior faces
3. **Compact face indexing**: Single array encodes all connectivity with signed integers
4. **Pre-sorted communication**: RemoteSendIndices ensure natural consumption order
5. **Minimal runtime logic**: Simple signed comparison determines all face handling

This design eliminates indirect addressing during the performance-critical flux computation phase, trading preprocessing complexity for runtime efficiency. The separation of connectivity from boundary conditions provides flexibility while maintaining performance.
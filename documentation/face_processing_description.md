# Face Buffer Runtime Operation and Boundary Processing

## Overview

The face buffer system provides a solution for managing face data in
partitioned meshes, enabling efficient computation of numerical fluxes through
Riemann solvers. The design achieves high performance by organizing data access
patterns to match the natural flow of computation, eliminating the need for
complex indexing during runtime.

## Face Buffer Structure and Natural Ordering

### The M Array

At the heart of the face buffer system is the M array, which stores all face
point data for a partition in a specific order called "natural traversal order."
This ordering follows the nested loop structure:

- Outermost loop: Elements (0 to K-1)
- Middle loop: Faces per element (0 to 3 for tetrahedra)
- Innermost loop: Points per face (depends on polynomial order)

This creates a linear sequence where face points are accessed in a predictable
pattern. For example, in a mesh with 100 elements, 4 faces per element, and 6
points per face, the M array contains 2,400 face points arranged sequentially.

### Face Classification

During the build phase, each position in the M array is classified as one of
three types:

1. **Boundary faces**: Face points on the domain boundary requiring boundary
   conditions
2. **Interior faces**: Face points connecting to other elements within the same
   partition
3. **Remote faces**: Face points connecting to elements in different partitions

This classification is stored in a parallel array that travels alongside the M
array.

## Runtime Traversal for Riemann Solvers

### The Fundamental Operation

When computing numerical fluxes, the Riemann solver requires two values at each
face point:

- The "minus" value (M): from the current element
- The "plus" value (P): from the neighboring element

The face buffer system retrieves these values efficiently during a single linear
traversal of the M array.

### Sequential Access Pattern

The runtime traversal follows this pattern:

1. Start at the beginning of the M array (position 0)
2. For each M position in sequence:
    - Retrieve the M value directly from the current position
    - Determine how to retrieve the corresponding P value based on face type
    - Compute the Riemann flux using M and P values
    - Advance to the next position

This sequential access pattern is cache-friendly and vectorization-friendly, as
it processes data in the order it's stored in memory.

## Retrieving P Values

The elegance of the face buffer design lies in how P values are retrieved for
different face types:

### Interior Faces (Local Neighbors)

For faces connecting to other elements within the same partition:

- A precomputed index array maps each interior M position to its corresponding P
  position
- The P value is retrieved from the same M array using this index
- Example: If M position 247 connects to position 1,891 in the same partition,
  the P value comes from M[1,891]

The system maintains a counter for interior faces that advances each time an
interior face is encountered, ensuring the correct index is used from the
precomputed array.

### Remote Faces (Inter-partition Connections)

For faces connecting to elements in other partitions:

- P values come from receive buffers populated by neighboring partitions
- Each partition sends its relevant M values to neighbors before the traversal
  begins
- The receiving partition stores these values in compact buffers, one per
  neighboring partition

The key insight is that remote P values are accessed in the exact order they're
needed during traversal:

1. Each remote partition has its own counter starting at 0
2. When a remote face from partition X is encountered, the next value is taken
   from partition X's receive buffer
3. The counter for partition X advances by one
4. This continues sequentially through the receive buffer

This design eliminates the need for index lookups or data reorganization—values
arrive and are consumed in natural traversal order.

### Boundary Faces

For domain boundary faces:

- The boundary condition type is stored alongside the face classification
- P values are computed based on the boundary condition formula
- Common examples include:
    - Wall boundaries: P = f(M) based on no-slip or slip conditions
    - Inflow boundaries: P = prescribed inflow values
    - Outflow boundaries: P = M (zero gradient)

## Communication Pattern for Remote Faces

### Pre-traversal Communication

Before the Riemann solver traversal begins:

1. Each partition identifies which of its M values are needed by neighboring
   partitions
2. These values are copied into send buffers in a specific order
3. Send buffers are exchanged between partitions (typically using MPI)
4. Each partition stores received values in partition-specific receive buffers

### Natural Order Preservation

The benefit of this approach is that both sending and receiving preserve
natural traversal order:

- Senders pack values in the order their neighbors will traverse them
- Receivers unpack values into buffers that will be consumed sequentially
- No reordering or indexing is needed during the actual computation

## Boundary Processing in Partitioned Meshes

### Unchanged from Single-Partition Meshes

Boundary condition processing in partitioned meshes remains identical to
unpartitioned meshes because:

1. **True boundaries are preserved**: Only actual domain boundaries are marked
   as boundary faces
2. **Partition interfaces are not boundaries**: Faces between partitions are
   classified as remote faces, not boundary faces
3. **BC application is local**: Each partition applies boundary conditions to
   its true boundary faces independently

### Boundary Condition Workflow

When the traversal encounters a boundary face:

1. The M value is retrieved from the current position
2. The boundary condition type is looked up from the classification array
3. The appropriate boundary condition formula computes P from M and possibly
   other data
4. The Riemann solver proceeds normally with these M and P values

This means existing boundary condition implementations work without modification
in the partitioned setting.

## Performance Benefits

This design achieves several performance advantages:

1. **Sequential memory access**: Both M array traversal and receive buffer
   access are purely sequential
2. **No indirect addressing during computation**: All indexing is resolved
   during the build phase
3. **Cache efficiency**: Data is accessed in the order it's stored
4. **Vectorization friendly**: Sequential access patterns enable SIMD operations
5. **Minimal branching**: Face type determines which code path, but within each
   type, operations are branch-free

## Example Traversal

Consider a partition processing element 15, face 2, point 3:

1. Compute M position: 15×4×6 + 2×6 + 3 = 375
2. Retrieve M value from position 375
3. Check face type at position 375:
    - If interior: Use interior counter (say 143) to get P from M array at
      precomputed position
    - If remote from partition 2: Use partition 2's counter (say 27) to get P
      from receive_buffer_2[27]
    - If boundary: Apply boundary condition formula to compute P
4. Compute Riemann flux with M and P
5. Advance all relevant counters
6. Continue to position 376

This process repeats for all 2,400 face points in perfect sequential order,
achieving optimal memory access patterns throughout the computation.

## Summary

The face buffer system transforms the complex problem of managing face data in
partitioned meshes into a simple linear traversal with sequential memory access.
By organizing data to match the natural flow of computation and pre-resolving
all indexing during the build phase, it is designed to achieve simplicity and
high performance. Boundary conditions integrate without special handling in 
the partitioned case, while inter-partition communication follows the same 
sequential access pattern as local operations.
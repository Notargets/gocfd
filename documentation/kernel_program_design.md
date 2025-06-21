# KernelProgram Design Document

## Overview

The KernelProgram system provides a high-level abstraction for generating and executing GPU/accelerator kernels for Discontinuous Galerkin (DG) methods using a **partition-parallel execution model**. The system divides large meshes into partitions that execute simultaneously on parallel hardware, with each partition processing its elements independently. A key innovation is the integration of the **M buffer** concept for efficient face data processing and inter-partition communication.

## Core Design Principles

1. **Partition-Parallel Execution**: Global mesh divided into Npart partitions that execute simultaneously on GPU blocks or CPU threads
2. **Variable Partition Sizes**: Each partition may have different element counts K[part]
3. **Pooled Memory Allocation**: Single allocation per array type with aligned partition offsets
4. **Tagged Alignment**: Different arrays use different alignment boundaries for performance
5. **OCCA-Native Parallelism**: Partitions map directly to OCCA @outer annotations
6. **M Buffer Architecture**: Face data organized in natural traversal order for sequential access
7. **Face-Based Communication**: Efficient inter-partition data exchange through face buffers

## Execution Model

### Mesh Partitioning

The global mesh with Kglobal elements is divided into Npart partitions with variable sizes:

```
Global Mesh (Kglobal elements)
    ├── Partition 0: K[0] elements
    ├── Partition 1: K[1] elements
    ├── ...
    └── Partition Npart-1: K[Npart-1] elements

where: Kglobal = sum(K[0] + K[1] + ... + K[Npart-1])
```

Partitioning is performed using METIS or similar graph partitioning libraries to ensure:
- **Load balance**: Each partition has approximately equal elements
- **Minimal edge cuts**: Reduced inter-partition communication
- **Partition count**: Chosen to match target hardware parallelism

### Parallel Execution

All partitions execute simultaneously on available hardware:
- **GPU**: Each partition runs on a separate CUDA block / OpenCL work-group
- **CPU**: Each partition runs on a separate OpenMP thread
- **Vectorization**: Within each partition, operations are vectorized over nodes

The partition count is configured to match the target hardware:
- **GPU**: Npart = 2-4× number of streaming multiprocessors (80-512)
- **CPU**: Npart = number of CPU threads (16-128)

## Face Buffer System

### M Buffer Concept

The M buffer stores all face point data for a partition in **natural traversal order**, following the nested loop structure:

```
for element in 0..K-1:
    for face in 0..NFACES-1:
        for point in 0..NFP-1:
            M[index++] = face_data[element][face][point]
```

This creates a linear array where face points are accessed sequentially, enabling:
- Cache-friendly memory access patterns
- Vectorization opportunities
- Minimal branching during traversal

### Face Classification

Each position in the M buffer is classified during the build phase:

```c
typedef enum {
    BOUNDARY_FACE = 0,  // Domain boundary requiring BC
    INTERIOR_FACE = 1,  // Connection within partition
    REMOTE_FACE = 2     // Connection to other partition
} FaceType;
```

### Local P Indexing

For interior faces connecting elements within the same partition, precomputed indices map M positions to their corresponding P positions:

```c
// For each interior M position, store where to find P in the same M buffer
int* localPIndices;  // Size: number of interior faces
int interiorCounter = 0;

// During traversal:
if (faceType[mIdx] == INTERIOR_FACE) {
    real_t M = faceBuffer[mIdx];
    real_t P = faceBuffer[localPIndices[interiorCounter++]];
    // Compute flux...
}
```

### Remote Buffer Exchange

For faces connecting to elements in other partitions, a carefully orchestrated communication pattern ensures sequential access without indexing:

**Key Principle**: The face buffer builder pre-computes all indexing. The ONLY memcpy is the MPI communication itself.

**Phase 1: Build-time Index Generation (face_buffer_builder.go)**
```go
// BuildFaceBuffer creates all indexing structures during mesh setup
type FaceBuffer struct {
// For interior faces within partition
LocalPIndices []uint32  // P position in M buffer for each interior point

// For remote faces to other partitions  
RemoteSendIndices map[uint32][]uint32 // [partitionID] -> which M indices to send
}

// The builder determines:
// 1. LocalPIndices: Where to find P values for interior faces
// 2. RemoteSendIndices: Which M values each neighbor needs as P values
// Both are computed to match the natural traversal order
```

**Phase 2: Direct MPI Communication**
```c
// THE ONLY MEMCPY: Direct MPI transfer from M buffer to receive buffers
// No intermediate buffer assembly needed - MPI reads directly from M using precomputed indices

// For each partition pair with connections:
for (int neighbor = 0; neighbor < Npart; ++neighbor) {
    if (hasConnection[part][neighbor]) {
        // MPI sends M values directly using RemoteSendIndices
        // These arrive in the receiver's natural traversal order
        MPI_Send_Indexed(M_buffer + M_offset[part],
                        remoteSendIndices[part][neighbor],
                        sendCount[part][neighbor],
                        destination: neighbor);
    }
}

// Receive buffers are populated in natural traversal order
MPI_Recv(recvBuffers[part][sender], recvCount[part][sender], sender);

**Phase 3: Sequential Traversal During Flux Computation**
```c
// Each partition maintains Npart-1 counters for sequential buffer traversal
int remoteCounters[NPART] = {0};  // One per potential neighbor

// During the natural M array traversal
for (int mIdx = 0; mIdx < totalFacePoints; ++mIdx) {
    real_t M = M_buffer[mIdx];
    real_t P;
    
    switch(faceTypes[mIdx]) {
        case INTERIOR_FACE:
            // Use precomputed local P index (from face_buffer_builder)
            P = M_buffer[localPIndices[interiorCounter++]];
            break;
            
        case REMOTE_FACE:
            // Which partition has this face's P value?
            int neighbor = faceNeighbor[mIdx];
            
            // Access the next P value from that partition's receive buffer
            // NO INDEXING - pure sequential access via counter
            P = recvBuffers[neighbor][remoteCounters[neighbor]++];
            break;
            
        case BOUNDARY_FACE:
            // Apply boundary condition
            P = applyBC(M, bcType[mIdx]);
            break;
    }
    
    // Compute flux with M and P
}
```
```

## OCCA Kernel Structure

Kernels use partition-level parallelism with OCCA-compliant loop nesting:

```c
@kernel void computeRHS(
    const int Npart,              // Number of partitions
    const int* K,                 // Elements per partition (variable)
    const real_t* U,              // Pooled solution array
    const int* U_offsets,         // Partition offsets into U
    const real_t* M_buffer,       // Face data in natural order
    const int* M_offsets,         // Partition offsets
    const FaceType* faceTypes,    // Face classifications
    const int* faceTypes_offsets,
    const int* localPIndices,     // Interior P positions
    const int* localP_offsets,
    const real_t* recvBuffers,    // Remote P values (post-exchange)
    const int* recv_offsets,
    real_t* RHS,                  // Pooled output array
    const int* RHS_offsets        // Partition offsets
) {
    // Each partition executes on its own GPU block / CPU thread
    for (int part = 0; part < Npart; ++part; @outer(0)) {
        
        // Get this partition's data pointers
        const real_t* myU = U + U_offsets[part];
        const real_t* myM = M_buffer + M_offsets[part];
        const FaceType* myFaceTypes = faceTypes + faceTypes_offsets[part];
        const int* myLocalP = localPIndices + localP_offsets[part];
        real_t* myRHS = RHS + RHS_offsets[part];
        const int myK = K[part];
        
        // Initialize per-partition counters for sequential buffer access
        int interiorCounter = 0;
        int remoteCounters[MAX_NEIGHBORS] = {0};  // One counter per neighbor partition
        
        // Process face contributions
        for (int elem = 0; elem < myK; ++elem) {
            for (int face = 0; face < NFACES; ++face) {
                
                // Vectorize over face points
                for (int fp = 0; fp < NFP; ++fp; @inner(0)) {
                    int mIdx = elem*NFACES*NFP + face*NFP + fp;
                    real_t M = myM[mIdx];
                    real_t P;
                    
                    // Retrieve P based on face type
                    switch(myFaceTypes[mIdx]) {
                        case BOUNDARY_FACE:
                            P = applyBC(M, bcType[mIdx]);
                            break;
                            
                        case INTERIOR_FACE:
                            // Use precomputed local P index from face_buffer_builder
                            P = myM[myLocalP[interiorCounter]];
                            interiorCounter++;
                            break;
                            
                        case REMOTE_FACE:
                            // Identify which partition has this face's P value
                            int neighbor = faceNeighbor[mIdx];
                            // Access next P value from that partition's receive buffer
                            // Pure sequential access - no indexing, just counter increment
                            // recvBuffers[part][neighbor] points to this partition's buffer from neighbor
                            P = recvBuffers[part][neighbor][remoteCounters[neighbor]];
                            remoteCounters[neighbor]++;  // Advance to next P value
                            break;
                    }
                    
                    // Compute numerical flux
                    real_t flux = riemannSolver(M, P, normal[mIdx]);
                    
                    // Accumulate to RHS
                    for (int n = 0; n < NP; ++n) {
                        myRHS[elem*NP + n] += LIFT[n][face*NFP + fp] * flux;
                    }
                }
            }
        }
    }
}
```

**Critical OCCA Constraints**:
1. `@inner` loops must be immediately nested within `@outer` loops
2. No intermediate loops allowed between `@outer` and `@inner`
3. All loops over elements must be inside the `@inner` loop
4. Multiple `@inner` loops are allowed but must have the same iteration count

## Memory Management

### Alignment Strategy

Different array types require different alignment boundaries for optimal performance:

```go
type AlignmentType int

const (
NoAlignment      AlignmentType = 1       // Natural alignment
CacheLineAlign   AlignmentType = 64      // Prevent false sharing between partitions
WarpAlign        AlignmentType = 128     // GPU warp-aligned access patterns
PageAlign        AlignmentType = 4096    // OS page boundary
LargePageAlign   AlignmentType = 2097152 // 2MB large pages for huge datasets
)
```

### Memory Specification

Each array specifies its size calculation and required alignment:

```go
type MemorySpec struct {
    Name      string
    SizeFunc  func(K, Np, Nfp, Nfaces int) int
    Alignment AlignmentType
}

var memoryLayout = []MemorySpec{
    // Solution arrays - cache line aligned to prevent false sharing
    {"U",          func(K, Np, _, _ int) int { return K * Np * 8 },           CacheLineAlign},
    {"RHS",        func(K, Np, _, _ int) int { return K * Np * 8 },           CacheLineAlign},
    
    // Face buffers - warp aligned for coalesced GPU access
    {"M_buffer",   func(K, _, Nfp, Nfaces int) int { return K * Nfaces * Nfp * 8 }, WarpAlign},
    {"faceTypes",  func(K, _, Nfp, Nfaces int) int { return K * Nfaces * Nfp * 1 }, NoAlignment},
    
    // Communication buffers - page aligned for efficient memcpy
    {"sendBuffer", func(K, _, Nfp, Nfaces int) int { return K * Nfaces * Nfp * 8 }, PageAlign},
    {"recvBuffer", func(K, _, Nfp, Nfaces int) int { return K * Nfaces * Nfp * 8 }, PageAlign},
    
    // Index arrays - no special alignment needed
    {"localPIndices", func(K, _, Nfp, Nfaces int) int { return K * Nfaces * Nfp * 4 }, NoAlignment},
    {"K",            func(_, _, _, _ int) int { return 8 }, NoAlignment},
}
```

### Pooled Allocation

Memory is allocated in large pools with proper alignment:

```go
func (kp *KernelProgram) AllocatePooledMemory() error {
    for _, spec := range memoryLayout {
        // Calculate total size across all partitions
        totalSize := 0
        offsets := make([]int, kp.NumPartitions)
        
        for part := 0; part < kp.NumPartitions; part++ {
            // Align offset to required boundary
            totalSize = alignUp(totalSize, spec.Alignment)
            offsets[part] = totalSize
            
            // Add this partition's data
            K := kp.K[part]
            partSize := spec.SizeFunc(K, kp.Np, kp.Nfp, kp.Nfaces)
            totalSize += partSize
        }
        
        // Allocate single contiguous buffer
        memory := kp.Device.Malloc(totalSize, spec.Alignment)
        kp.pooledMemory[spec.Name] = memory
        
        // Store offsets for kernel use
        offsetsMem := kp.Device.Malloc(kp.NumPartitions * 8)
        offsetsMem.CopyFrom(offsets)
        kp.pooledMemory[spec.Name + "_offsets"] = offsetsMem
    }
    
    return nil
}
```

### Memory Layout Example

For 3 partitions with K=[100, 150, 120] elements and Np=20, Nfp=6:

```
M_buffer (warp-aligned at 128 bytes):
[Partition 0: 100*4*6*8 = 19200 bytes][pad 64 bytes][Partition 1: 150*4*6*8 = 28800 bytes][pad 0 bytes][Partition 2: 120*4*6*8 = 23040 bytes]
 ^                                                    ^                                                      ^
 offset[0] = 0                                        offset[1] = 19264                                    offset[2] = 48064

recvBuffer (page aligned at 4096 bytes):
[Partition 0: buffer][pad to 4KB][Partition 1: buffer][pad to 4KB][Partition 2: buffer]
 ^                                ^                                 ^
 offset[0] = 0                    offset[1] = 20480                offset[2] = 49152
```

## Code Generation

### Preamble Generation

The kernel preamble includes essential definitions:

```c
// Generated preamble
#define ORDER @Order@             // Polynomial order
#define NP @Np@                   // Nodes per element
#define NFP @Nfp@                 // Face nodes per element  
#define NFACES 4                  // Faces per tetrahedral element
#define NPART @NumPartitions@     // Total number of partitions

// Type definitions
typedef @FloatType@ real_t;       // double or float
typedef @IntType@ int_t;          // int64_t or int32_t

// Face types
typedef enum {
    BOUNDARY_FACE = 0,
    INTERIOR_FACE = 1,
    REMOTE_FACE = 2
} FaceType;

// Static matrices (embedded as constants)
const real_t Dr[@Np@][@Np@] = { ... };
const real_t Ds[@Np@][@Np@] = { ... };
const real_t Dt[@Np@][@Np@] = { ... };
const real_t LIFT[@Np@][@Nfp@*@Nfaces@] = { ... };

// Note: No K definition - it varies per partition
```

## Kernel Execution

### RunKernel Implementation

The kernel execution relies on the kernel's own parallelism specification:

```go
func (kp *KernelProgram) RunKernel(name string, args ...interface{}) error {
    kernel, exists := kp.kernels[name]
    if !exists {
        return fmt.Errorf("kernel %s not found", name)
    }
    
    // Let the kernel control its own parallelism
    // The @outer loop in the kernel specifies partition parallelism
    
    return kernel.RunWithArgs(args...)
}
```

## Design Benefits

1. **Variable Partition Support**: Natural handling of non-uniform mesh partitions through K array
2. **Single Allocation Per Array**: Avoids memory fragmentation, guarantees contiguity
3. **Optimal Alignment**: Each array type gets appropriate alignment for its access pattern
4. **OCCA Compliance**: Proper loop nesting without conflicting parallelism specifications
5. **Clean Kernel Interface**: Simple base_pointer + offset[part] indexing
6. **Sequential Memory Access**: M buffer design enables cache-friendly traversal patterns
7. **Minimal Data Movement**: The ONLY operation is indexed memcpy within shared memory - no intermediate buffers needed
8. **Pre-computed Indexing**: face_buffer_builder.go generates all indices at build time:
   - LocalPIndices for interior face P lookups
   - RemoteSendIndices for copying M values between partitions
9. **Index-Free Runtime Access**: Received P values are consumed in natural traversal order using simple counters
10. **Performance Optimizations**:
   - Cache-line alignment prevents false sharing between partitions
   - Memory alignment optimizes coalesced access for memcpy operations
   - Warp alignment optimizes GPU coalesced access patterns
   - Sequential access patterns enable vectorization
   - No indirect memory access during flux computation
   - Counter-based traversal eliminates all runtime index lookups
   - Shared memory access avoids any network or bus communication

## Usage Example

```go
// Create kernel program with variable partitions
partitionElementCounts := []int{102, 98, 105, 97, 103, ...} // From METIS
kp := NewKernelProgram(device, Config{
NumPartitions: 64,
K:             partitionElementCounts,
Order:         4,
FloatType:     Float64,
})

// Generate static data and type definitions
kp.AddStaticMatrix("Dr", Dr)
kp.AddStaticMatrix("Ds", Ds)
kp.AddStaticMatrix("Dt", Dt)
kp.GenerateKernelPreamble()

// Build kernels
kp.BuildKernel(rhsKernelSource, "computeRHS")
kp.BuildKernel(timeStepKernel, "updateSolution")

// Allocate pooled memory with proper alignment
kp.AllocatePooledMemory()

// Build face buffers for each partition using face_buffer_builder.go
faceBuffers := make([]*FaceBuffer, nPart)
for part := 0; part < nPart; part++ {
fb, err := BuildFaceBuffer(elements[part])
if err != nil {
log.Fatal(err)
}
faceBuffers[part] = fb

// Upload precomputed indices to device
kp.SetMemory(fmt.Sprintf("localPIndices_%d", part), fb.LocalPIndices)
kp.SetMemory(fmt.Sprintf("remoteSendIndices_%d", part), fb.RemoteSendIndices)
}

// Time stepping loop
for step := 0; step < nSteps; step++ {
// THE ONLY MEMCPY: Copy M values between partitions using precomputed indices
// This kernel performs indexed memcpy from M buffers to receive buffers
// within shared memory - no MPI or network communication
kp.RunKernel("exchangeFaceData",
kp.NumPartitions,
kp.GetMemory("M_buffer"), kp.GetMemory("M_offsets"),
kp.GetMemory("recvBuffer"), kp.GetMemory("recvBuffer_offsets"),
kp.GetMemory("remoteSendIndices"),
kp.GetMemory("sendCounts"))

// Compute RHS with purely sequential traversal
// Interior P values use LocalPIndices (from face_buffer_builder)
// Remote P values are accessed from receive buffers using counters (no indices)
kp.RunKernel("computeRHS",
kp.NumPartitions,
kp.GetMemory("K"),
kp.GetMemory("U"), kp.GetMemory("U_offsets"),
kp.GetMemory("M_buffer"), kp.GetMemory("M_offsets"),
kp.GetMemory("faceTypes"), kp.GetMemory("faceTypes_offsets"),
kp.GetMemory("localPIndices"), kp.GetMemory("localP_offsets"),
kp.GetMemory("recvBuffer"), kp.GetMemory("recvBuffer_offsets"),
kp.GetMemory("RHS"), kp.GetMemory("RHS_offsets"))

// Time integration...
}
```

This design enables high-performance DG computations while maintaining code clarity and portability across different parallel backends. The M buffer architecture ensures optimal memory access patterns, while the face classification system enables efficient handling of boundaries and inter-partition communication.
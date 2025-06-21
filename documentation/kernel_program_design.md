# KernelProgram Design Document

## Overview

The KernelProgram system provides a high-level abstraction for generating and executing GPU/accelerator kernels for Discontinuous Galerkin (DG) methods using a **partition-parallel execution model**. The system divides large meshes into partitions that execute simultaneously on parallel hardware within a shared memory architecture. A key innovation is the integration of the **M buffer** concept for efficient face data processing and inter-partition communication.

## Core Design Principles

1. **Partition-Parallel Execution**: Global mesh divided into Npart partitions that execute simultaneously on GPU blocks or CPU threads
2. **Variable Partition Sizes**: Each partition may have different element counts; kArray[part] gives elements for partition part
3. **Pooled Memory Allocation**: Single allocation per array type with aligned partition offsets
4. **Tagged Alignment**: Different arrays use different alignment boundaries for performance
5. **OCCA-Native Parallelism**: Partitions map directly to OCCA @outer annotations
6. **M Buffer Architecture**: Face data organized in natural traversal order for sequential access
7. **Face-Based Communication**: Efficient inter-partition data exchange through face buffers in shared memory

## Execution Model

### Mesh Partitioning

The global mesh with Kglobal elements is divided into Npart partitions with variable sizes:

```
Global Mesh (Kglobal elements)
    ├── Partition 0: kArray[0] elements
    ├── Partition 1: kArray[1] elements
    ├── ...
    └── Partition Npart-1: kArray[Npart-1] elements

where: Kglobal = sum(kArray[0] + kArray[1] + ... + kArray[Npart-1])
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
for element in 0..K-1:    // K is the local element count for this partition
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

### Face Buffer Builder Integration

The face_buffer_builder.go creates essential indexing structures:

```go
type FaceBuffer struct {
    // Dimensions
    Nfp             uint32 // Face points per face
    Nfaces          uint32 // Faces per element (4 for tetrahedra)
    K               uint32 // Number of elements in this partition (local K value)
    TotalFacePoints uint32 // Total M buffer size = K * Nfaces * Nfp

    // Face classification
    FaceTypes []FaceType // Type of each face point connection

    // Local interior P indices
    LocalPIndices []uint32 // P position in M buffer for each interior point

    // Remote partition send indices
    RemoteSendIndices map[uint32][]uint32 // [partitionID] -> indices into M buffer
}
```

These structures enable:
- **LocalPIndices**: Direct lookup of P values for interior faces without runtime computation
- **RemoteSendIndices**: Identifies which M values to copy for neighboring partitions

## OCCA Memory Management

### Memory Allocation

The system uses OCCA's device memory allocation API:

```go
func (kp *KernelProgram) AllocatePooledMemory() error {
    // First, allocate the kArray (element counts per partition)
    // This is a global array, not per-partition
    kArrayData := make([]int32, kp.NumPartitions)
    for i := 0; i < kp.NumPartitions; i++ {
        kArrayData[i] = int32(kp.ElementsPerPartition[i])
    }
    kp.memory["kArray"] = kp.device.Malloc(
        int64(kp.NumPartitions * 4), // 4 bytes per int32
        unsafe.Pointer(&kArrayData[0]),
        nil,
    )
    
    // Calculate total sizes for each per-partition array type
    for _, spec := range memoryLayout {
        totalSize := 0
        offsets := make([]int, kp.NumPartitions)
        
        for part := 0; part < kp.NumPartitions; part++ {
            // Align offset to required boundary
            totalSize = alignUp(totalSize, spec.Alignment)
            offsets[part] = totalSize
            
            // Add this partition's data
            K := kp.ElementsPerPartition[part]  // Local K value for this partition
            partSize := spec.SizeFunc(K, kp.Np, kp.Nfp, kp.Nfaces)
            totalSize += partSize
        }
        
        // Allocate using OCCA API
        memory := kp.device.Malloc(int64(totalSize), nil, nil)
        kp.memory[spec.Name] = memory
        
        // Store offsets in device memory
        offsetsMem := kp.device.Malloc(int64(kp.NumPartitions * 8), nil, nil)
        offsetsMem.CopyFrom(unsafe.Pointer(&offsets[0]), int64(len(offsets)*8))
        kp.memory[spec.Name + "_offsets"] = offsetsMem
    }
    
    return nil
}
```

### Alignment Strategy

Different array types require different alignment boundaries:

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

```go
var memoryLayout = []MemorySpec{
    // Solution arrays - cache line aligned to prevent false sharing
    // Note: K parameter in SizeFunc is the partition-local element count
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
}
```

## OCCA Kernel Structure

### OKL Kernel with Explicit Parallelism

Kernels use OCCA's @outer/@inner parallelism model:

```c
@kernel void computeRHS(
    const int Npart,              // Number of partitions
    const int* kArray,            // Array of elements per partition (kArray[0], kArray[1], ..., kArray[Npart-1])
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
        const int K = kArray[part];  // Extract partition-local element count
        
        // Initialize per-partition counters for sequential buffer access
        int interiorCounter = 0;
        int remoteCounters[MAX_NEIGHBORS] = {0};  // One counter per neighbor partition
        
        // Process face contributions
        for (int elem = 0; elem < K; ++elem) {
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

## Kernel Execution

### Building Kernels

```go
func (kp *KernelProgram) BuildKernel(kernelSource, kernelName string) (*gocca.OCCAKernel, error) {
    // Combine preamble with kernel source
    fullSource := kp.kernelPreamble + "\n" + kernelSource
    
    // Build kernel using OCCA API
    kernel, err := kp.device.BuildKernelFromString(fullSource, kernelName, nil)
    if err != nil {
        return nil, fmt.Errorf("failed to build kernel %s: %w", kernelName, err)
    }
    
    kp.RegisterKernel(kernelName, kernel)
    return kernel, nil
}
```

### Running Kernels

For OKL kernels, the parallelism is specified in the kernel source:

```go
func (kp *KernelProgram) RunKernel(name string, args ...interface{}) error {
    kernel, exists := kp.kernels[name]
    if !exists {
        return fmt.Errorf("kernel %s not found", name)
    }
    
    // For OKL kernels, parallelism is controlled by @outer/@inner in kernel source
    // For non-OKL kernels, use SetRunDims as shown in existing implementation
    
    return kernel.RunWithArgs(args...)
}
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
#define MAX_NEIGHBORS 64          // Maximum possible neighbor partitions

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
__constant__ real_t Dr[@Np@][@Np@] = { ... };
__constant__ real_t Ds[@Np@][@Np@] = { ... };
__constant__ real_t Dt[@Np@][@Np@] = { ... };
__constant__ real_t LIFT[@Np@][@Nfp@*@Nfaces@] = { ... };

// Note: K (elements per partition) is passed as the runtime array kArray since each partition has different element counts
```

## Communication Architecture

The face buffer exchange system uses a two-phase approach for optimal performance:

1. **Build-Time Index Generation**: face_buffer_builder.go creates:
   - **LocalPIndices**: Maps each interior M position to its P position within the partition
   - **RemoteSendIndices**: For each neighbor partition, lists which M values to gather and send

2. **Runtime Data Movement** (two kernels):
   - **gatherSendBuffers**: Each partition gathers its M values into contiguous send buffers using RemoteSendIndices
   - **copySendToReceiveBuffers**: Send buffers are copied to receive buffers in shared memory

3. **Natural Order Preservation**: RemoteSendIndices are ordered to match the receiving partition's traversal order, enabling sequential access of receive buffers.

4. **Sequential Runtime Access**: During flux computation:
   - Interior P values: Look up using LocalPIndices[interiorCounter++]
   - Remote P values: Access using recvBuffers[neighbor][remoteCounter[neighbor]++]

This design creates contiguous buffers for efficient memcpy and eliminates all runtime indexing during the performance-critical flux computation phase.

**Note**: In a shared memory system, we could potentially skip the send-to-receive buffer copy by having partitions read directly from each other's send buffers. This optimization is left for future implementation.

## Usage Example

```go
type Config struct {
    Order                int
    NumPartitions        int   // Number of partitions (Npart)
    ElementsPerPartition []int // Array of element counts per partition (kArray)
    FloatType           DataType
    IntType             DataType
}

// Generate static data and type definitions
kp.AddStaticMatrix("Dr", Dr)
kp.AddStaticMatrix("Ds", Ds)
kp.AddStaticMatrix("Dt", Dt)
kp.GenerateKernelPreamble()

// Build kernels
kp.BuildKernel(gatherKernel, "gatherSendBuffers")
kp.BuildKernel(copyKernel, "copySendToReceiveBuffers") 
kp.BuildKernel(rhsKernelSource, "computeRHS")

// Allocate pooled memory with proper alignment
kp.AllocatePooledMemory()

// Build face buffers for each partition
faceBuffers := make([]*FaceBuffer, nPart)
localPIndexArrays := make([]*gocca.OCCAMemory, nPart)
remoteSendArrays := make(map[uint32][]*gocca.OCCAMemory)

for part := 0; part < nPart; part++ {
    fb, err := BuildFaceBuffer(elements[part])
    if err != nil {
        log.Fatal(err)
    }
    faceBuffers[part] = fb
    
    // Upload LocalPIndices to device
    localPMem := device.Malloc(int64(len(fb.LocalPIndices)*4), 
                               unsafe.Pointer(&fb.LocalPIndices[0]), nil)
    localPIndexArrays[part] = localPMem
    
    // Upload RemoteSendIndices for each neighbor
    for neighborID, indices := range fb.RemoteSendIndices {
        sendMem := device.Malloc(int64(len(indices)*4),
                                unsafe.Pointer(&indices[0]), nil)
        remoteSendArrays[neighborID] = append(remoteSendArrays[neighborID], sendMem)
    }
    
    // Upload face types
    faceTypesMem := device.Malloc(int64(len(fb.FaceTypes)), 
                                 unsafe.Pointer(&fb.FaceTypes[0]), nil)
    kp.SetMemory(fmt.Sprintf("faceTypes_%d", part), faceTypesMem)
}

// Time stepping loop
for step := 0; step < nSteps; step++ {
    // Phase 1: Copy M values between partitions using precomputed indices
    kp.RunKernel("exchangeFaceData",
        kp.NumPartitions,
        kp.GetMemory("M_buffer"), kp.GetMemory("M_offsets"),
        kp.GetMemory("recvBuffer"), kp.GetMemory("recvBuffer_offsets"),
        remoteSendArrays,
        sendCounts)
    
    // Phase 2: Compute RHS with purely sequential traversal
    kp.RunKernel("computeRHS",
        kp.NumPartitions,
        kp.GetMemory("kArray"),  // Array of elements per partition
        kp.GetMemory("U"), kp.GetMemory("U_offsets"),
        kp.GetMemory("M_buffer"), kp.GetMemory("M_offsets"),
        kp.GetMemory("faceTypes"), kp.GetMemory("faceTypes_offsets"),
        localPIndexArrays, localPOffsets,
        kp.GetMemory("recvBuffer"), kp.GetMemory("recvBuffer_offsets"),
        kp.GetMemory("RHS"), kp.GetMemory("RHS_offsets"))
    
    // Time integration...
}
```

## Design Benefits

1. **OCCA Compliance**: Uses standard OCCA memory allocation and kernel execution APIs
2. **Variable Partition Support**: Natural handling of non-uniform mesh partitions through kArray
3. **Single Allocation Per Array**: Avoids memory fragmentation, uses OCCA's pooled allocation
4. **Optimal Alignment**: Each array type gets appropriate alignment for its access pattern
5. **OKL Parallelism**: Proper @outer/@inner loop nesting for GPU/CPU portability
6. **Sequential Memory Access**: M buffer design enables cache-friendly traversal patterns
7. **Minimal Data Movement**: Only indexed memcpy within shared memory - no network overhead
8. **Pre-computed Indexing**: face_buffer_builder.go generates all indices at build time
9. **Index-Free Runtime Access**: Received P values consumed in natural order using counters
10. **Performance Optimizations**:
   - Cache-line alignment prevents false sharing between partitions
   - Memory alignment optimizes coalesced access for memcpy operations
   - Warp alignment optimizes GPU coalesced access patterns
   - Sequential access patterns enable vectorization
   - No indirect memory access during flux computation
   - Shared memory architecture eliminates communication overhead

This design enables high-performance DG computations while maintaining code clarity and full compatibility with the OCCA framework.
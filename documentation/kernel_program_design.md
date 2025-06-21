# KernelProgram Design Document

## Overview

The KernelProgram system provides a high-level abstraction for generating and executing GPU/accelerator kernels for Discontinuous Galerkin (DG) methods using a **partition-parallel execution model**. The system divides large meshes into partitions that execute simultaneously on parallel hardware, with each partition processing its elements independently.

## Core Design Principles

1. **Partition-Parallel Execution**: Global mesh divided into Npart partitions that execute simultaneously on GPU blocks or CPU threads
2. **Variable Partition Sizes**: Each partition may have different element counts K[part]
3. **Pooled Memory Allocation**: Single allocation per array type with aligned partition offsets
4. **Tagged Alignment**: Different arrays use different alignment boundaries for performance
5. **OCCA-Native Parallelism**: Partitions map directly to OCCA @outer annotations without external dimension control

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

### OCCA Kernel Structure

Kernels use partition-level parallelism with OCCA-compliant loop nesting:

```c
@kernel void computeRHS(
    const int Npart,              // Number of partitions
    const int* K,                 // Elements per partition (variable)
    const real_t* U,              // Pooled input array
    const int* U_offsets,         // Partition offsets into U
    const real_t* geoFactors,     // Pooled geometric factors
    const int* geo_offsets,       // Partition offsets
    real_t* RHS,                  // Pooled output array
    const int* RHS_offsets        // Partition offsets
) {
    // Each partition executes on its own GPU block / CPU thread
    for (int part = 0; part < Npart; ++part; @outer(0)) {
        
        // Get this partition's data pointers
        const real_t* myU = U + U_offsets[part];
        const real_t* myGeo = geoFactors + geo_offsets[part];
        real_t* myRHS = RHS + RHS_offsets[part];
        const int myK = K[part];
        
        // Vectorize over nodes within element
        // NOTE: @inner must be immediately nested within @outer
        for (int node = 0; node < NP; ++node; @inner(0)) {
            
            // Process all elements in this partition
            // Element loop MUST be inside @inner
            for (int elem = 0; elem < myK; ++elem) {
                
                // Partition-local index
                int idx = elem*NP + node;
                
                // Computation using partition-local arrays
                myRHS[idx] = computeLocalRHS(myU, myGeo, elem, node);
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
    DMAAlign         AlignmentType = 131072  // 128KB for parallel DMA transfers
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
    
    // Geometric factors - cache line aligned
    {"geoFactors", func(K, Np, _, _ int) int { return K * Np * 9 * 8 },       CacheLineAlign},
    
    // Face data - warp aligned for coalesced GPU access
    {"faceM",      func(K, _, Nfp, Nfaces int) int { return K * Nfaces * Nfp * 8 }, WarpAlign},
    {"faceP",      func(K, _, Nfp, Nfaces int) int { return K * Nfaces * Nfp * 8 }, WarpAlign},
    
    // Communication buffers - DMA aligned for fast parallel transfers
    {"sendBuffer", func(K, _, Nfp, Nfaces int) int { return K * Nfaces * Nfp * 8 }, DMAAlign},
    {"recvBuffer", func(K, _, Nfp, Nfaces int) int { return K * Nfaces * Nfp * 8 }, DMAAlign},
    
    // Integer arrays - no special alignment needed
    {"faceTypes",  func(K, _, _, Nfaces int) int { return K * Nfaces * 4 },         NoAlignment},
    {"K",          func(_, _, _, _ int) int { return 8 },                           NoAlignment}, // Per partition
}
```

### Pooled Allocation with Aligned Offsets

Memory allocation computes aligned offsets for each partition within a single allocation:

```go
func (kp *KernelProgram) AllocatePooledMemory() error {
    // Allocate K array (element counts per partition)
    kp.memory["K"] = kp.device.Malloc(int64(kp.NumPartitions * 8), 
                                      unsafe.Pointer(&kp.K[0]), nil)
    
    // Allocate pooled arrays with alignment
    for _, spec := range memoryLayout {
        if spec.Name == "K" {
            continue // Already allocated
        }
        
        offsets := make([]int64, kp.NumPartitions)
        currentOffset := int64(0)
        
        // Calculate aligned offset for each partition
        for p := 0; p < kp.NumPartitions; p++ {
            // Align current offset to specified boundary
            if spec.Alignment > 1 {
                alignment := int64(spec.Alignment)
                currentOffset = ((currentOffset + alignment - 1) / alignment) * alignment
            }
            
            offsets[p] = currentOffset
            
            // Calculate this partition's size
            size := spec.SizeFunc(kp.K[p], kp.Np, kp.Nfp, kp.Nfaces)
            currentOffset += int64(size)
        }
        
        // Total allocation = sum of all partition sizes + alignment padding
        totalSize := currentOffset
        
        // Single allocation for all partitions
        kp.memory[spec.Name] = kp.device.Malloc(totalSize, nil, nil)
        
        // Store offsets for kernel use (convert to appropriate type)
        offsetsArray := make([]int64, kp.NumPartitions)
        for i := range offsets {
            offsetsArray[i] = offsets[i] / kp.GetElementSize(spec.Name)
        }
        
        // Copy offsets to device
        offsetsMem := kp.device.Malloc(int64(kp.NumPartitions * 8), 
                                       unsafe.Pointer(&offsetsArray[0]), nil)
        kp.memory[spec.Name + "_offsets"] = offsetsMem
    }
    
    return nil
}
```

### Memory Layout Example

For 3 partitions with K=[100, 150, 120] elements and Np=20:

```
Array U (cache-line aligned at 64 bytes):
[Partition 0: 100*20*8 = 16000 bytes][pad 0 bytes][Partition 1: 150*20*8 = 24000 bytes][pad 0 bytes][Partition 2: 120*20*8 = 19200 bytes]
 ^                                                 ^                                                   ^
 offset[0] = 0                                     offset[1] = 16000                                  offset[2] = 40000

sendBuffer (DMA aligned at 128KB = 131072 bytes):
[Partition 0: 100*4*6*8 = 19200 bytes][pad 111872 bytes][Partition 1: 150*4*6*8 = 28800 bytes][pad 102272 bytes][Partition 2: ...]
 ^                                                       ^                                                        ^
 offset[0] = 0                                           offset[1] = 131072                                      offset[2] = 262144
```

## Code Generation

### Preamble Generation

The kernel preamble includes only essential definitions without assuming uniform K:

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
    
    // NO SetRunDims - let the kernel control its own parallelism
    // The @outer loop in the kernel specifies partition parallelism
    
    return kernel.RunWithArgs(args...)
}
```

## Face Data Exchange

The all-to-all face exchange pattern for inter-partition communication:

```c
@kernel void exchangeFaceData(
    const int Npart,
    const int* K,
    const real_t* sendBuffer,      // DMA-aligned for efficiency
    const int* send_offsets,
    real_t* recvBuffer,           // DMA-aligned for efficiency  
    const int* recv_offsets,
    const int* faceConnectivity,  // [part][face] -> (neighbor_part, neighbor_face)
    const int* faceCounts         // Number of faces per partition
) {
    // Each partition handles its receive operations
    for (int recvPart = 0; recvPart < Npart; ++recvPart; @outer(0)) {
        
        const int* myConn = faceConnectivity + recvPart * MAX_FACES;
        real_t* myRecv = recvBuffer + recv_offsets[recvPart];
        
        // Parallel copy across face points
        for (int fp = 0; fp < NFP; ++fp; @inner(0)) {
            
            // Process all faces for this partition
            for (int face = 0; face < faceCounts[recvPart]; ++face) {
                int sendPart = myConn[face * 2];
                int sendFace = myConn[face * 2 + 1];
                
                if (sendPart >= 0) {  // Valid connection
                    const real_t* theirSend = sendBuffer + send_offsets[sendPart];
                    
                    // Copy face data with coalesced access
                    for (int node = fp; node < NFP; node += NFP) {
                        int recvIdx = face * NFP + node;
                        int sendIdx = sendFace * NFP + node;
                        myRecv[recvIdx] = theirSend[sendIdx];
                    }
                }
            }
        }
    }
}
```

## Design Benefits

1. **Variable Partition Support**: Natural handling of non-uniform mesh partitions through K array
2. **Single Allocation Per Array**: Avoids memory fragmentation, guarantees contiguity
3. **Optimal Alignment**: Each array type gets appropriate alignment for its access pattern
4. **OCCA Compliance**: Proper loop nesting without conflicting parallelism specifications
5. **Clean Kernel Interface**: Simple base_pointer + offset[part] indexing
6. **Performance Optimizations**:
    - Cache-line alignment prevents false sharing between partitions
    - DMA alignment enables fast parallel memcpy for communication
    - Warp alignment optimizes GPU coalesced access patterns
    - No wasted memory on maximum-based allocation

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
kp.BuildKernel(exchangeKernelSource, "exchangeFaceData")

// Allocate pooled memory with proper alignment
kp.AllocatePooledMemory()

// Time stepping loop
for step := 0; step < nSteps; step++ {
    // Exchange face data between partitions
    kp.RunKernel("exchangeFaceData",
        kp.NumPartitions,
        kp.GetMemory("K"),
        kp.GetMemory("sendBuffer"), kp.GetMemory("sendBuffer_offsets"),
        kp.GetMemory("recvBuffer"), kp.GetMemory("recvBuffer_offsets"),
        kp.GetMemory("faceConnectivity"),
        kp.GetMemory("faceCounts"))
    
    // Compute RHS with all partitions in parallel
    kp.RunKernel("computeRHS",
        kp.NumPartitions,
        kp.GetMemory("K"),
        kp.GetMemory("U"), kp.GetMemory("U_offsets"),
        kp.GetMemory("geoFactors"), kp.GetMemory("geoFactors_offsets"),
        kp.GetMemory("RHS"), kp.GetMemory("RHS_offsets"))
    
    // Time integration...
}
```

This design enables high-performance DG computations while maintaining code clarity and portability across different parallel backends.
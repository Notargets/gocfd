# KernelProgram Design Document

## Overview

The KernelProgram system provides a high-level abstraction for generating and executing GPU/accelerator kernels for Discontinuous Galerkin (DG) methods using a **partition-parallel execution model**. The system divides large meshes into partitions that execute simultaneously on parallel hardware, with each partition processing its elements independently.

## Core Design Principles

1. **Partition-Parallel Execution**: Global mesh divided into Npart partitions that execute simultaneously on GPU blocks or CPU threads
2. **Element-Blocked Data Layout**: Within each partition, nodes of an element are contiguous for cache efficiency
3. **Partition-Blocked Memory**: Partitions are contiguous in memory to enable coalesced access across parallel units
4. **Static Data Embedding**: Small matrices (Dr, Ds, Dt, LIFT) compiled into kernels as constants
5. **OCCA-Native Parallelism**: Partitions map directly to OCCA @outer annotations for hardware parallelism

## Execution Model

### Mesh Partitioning

The global mesh with Kglobal elements is divided into Npart partitions:

```
Global Mesh (Kglobal elements)
    ├── Partition 0: K elements
    ├── Partition 1: K elements
    ├── ...
    └── Partition Npart-1: K elements

where: Kglobal = Npart × K (approximately)
```

### Parallel Execution

All partitions execute simultaneously:
- **GPU**: Each partition runs on a separate CUDA block / OpenCL work-group
- **CPU**: Each partition runs on a separate thread
- **Vectorization**: Within each partition, operations are vectorized over elements

### OCCA Kernel Structure

Kernels use partition-level parallelism as the primary @outer loop:

```c
@kernel void computeRHS(
    const int Npart,          // Number of partitions
    const int K,              // Elements per partition  
    const real_t* U,          // Input: [Npart][K][NP]
    const real_t* geoFactors, // Geometric factors: [Npart][K][NP]
    real_t* RHS               // Output: [Npart][K][NP]
) {
    // Each partition executes on its own GPU block / CPU thread
    for (int part = 0; part < Npart; ++part; @outer(0)) {
        
        // Process all elements in this partition
        for (int elem = 0; elem < K; ++elem) {
            
            // Vectorize over nodes within element
            for (int node = 0; node < NP; ++node; @inner(0)) {
                
                // Global index for partition-blocked data
                int idx = part*K*NP + elem*NP + node;
                
                // Partition-local operations
                // All data access is within this partition's memory region
            }
        }
    }
}
```

## Architecture

### KernelProgram Structure

```go
type KernelProgram struct {
    // Partition configuration
    NumPartitions   int      // Npart - number of partitions
    ElementsPerPart int      // K - elements per partition
    TotalElements   int      // Kglobal - total elements across all partitions
    
    // Element configuration
    Order           int      // Polynomial order (N)
    Np              int      // Nodes per element
    Nfp             int      // Face points per element face
    Nfaces          int      // Faces per element (4 for tet)
    
    // Data precision
    FloatType       DataType // Float32 or Float64
    IntType         DataType // Int32 or Int64
    
    // Static data (shared across all partitions)
    StaticMatrices  map[string]utils.Matrix
    
    // Generated code
    kernelPreamble  string
    
    // Runtime resources
    device          *gocca.OCCADevice
    kernels         map[string]*gocca.OCCAKernel
    memory          map[string]*gocca.OCCAMemory  // Partition-blocked arrays
}
```

### Memory Layout

#### Partition-Blocked Organization

All dynamic data uses partition-blocked layout:

```
Array U[Npart][K][NP] stored contiguously as:
[partition0_elem0_node0, ..., partition0_elem0_nodeNP-1,  // First element of first partition
 partition0_elem1_node0, ..., partition0_elem1_nodeNP-1,  // Second element
 ...
 partition0_elemK-1_node0, ...,                           // Last element of first partition
 partition1_elem0_node0, ...,                             // First element of second partition
 ...]

Index calculation: idx = part*K*NP + elem*NP + node
```

This layout ensures:
- **Coalesced access**: Adjacent partitions (GPU blocks) access adjacent memory regions
- **Cache efficiency**: Element nodes remain contiguous for element-local operations
- **Minimal conflicts**: Partitions work on separate memory regions

#### Memory Allocation Sizes

```go
nodeCount := Npart * K * Np           // Total nodes across all partitions
faceCount := Npart * K * Nfaces * Nfp // Total face nodes
```

## Code Generation

### Partition-Aware Preamble

The generated preamble includes partition-aware utilities:

```c
// Constants
#define NPART @NumPartitions@    // Total partitions
#define K @ElementsPerPart@       // Elements per partition
#define NP @Np@                   // Nodes per element
#define NFP @Nfp@                 // Face nodes per element
#define NFACES 4                  // Faces per tet

// Indexing macros for partition-blocked data
#define PART_OFFSET(part) ((part) * K * NP)
#define ELEM_OFFSET(part, elem) (PART_OFFSET(part) + (elem) * NP)
#define NODE_INDEX(part, elem, node) (ELEM_OFFSET(part, elem) + (node))

// Type definitions
typedef @FloatType@ real_t;
typedef @IntType@ int_t;
```

### Partition-Aware Matrix Operations

Generated matrix operations handle partition offsets:

```c
// Matrix multiplication for partition-blocked data
inline void matMul_Dr_Partition(
    const real_t* U,      // [Npart][K][NP]
    real_t* result,       // [Npart][K][NP]
    int part,             // Current partition ID
    int K                 // Elements in partition
) {
    int partOffset = PART_OFFSET(part);
    
    // Process all elements in this partition
    for (int elem = 0; elem < K; ++elem) {
        // Apply Dr to each element
        for (int i = 0; i < NP; ++i) {
            real_t sum = REAL_ZERO;
            #pragma unroll
            for (int j = 0; j < NP; ++j) {
                sum += Dr[i][j] * U[partOffset + elem*NP + j];
            }
            result[partOffset + elem*NP + i] = sum;
        }
    }
}

// Convenience function that extracts partition ID from kernel context
inline void matMul_Dr_Auto(const real_t* U, real_t* result, int K) {
    // In OCCA kernels, can get partition ID from loop index
    extern int part;  // Set by @outer loop
    matMul_Dr_Partition(U, result, part, K);
}
```

### Physical Derivative Operations

```c
// Compute derivatives for one partition
inline void computePhysDerivPartition(
    const real_t* U,          // Solution data
    const real_t* geoFactors, // rx,ry,rz,sx,sy,sz,tx,ty,tz
    real_t* DX,               // Output derivatives
    real_t* DY,
    real_t* DZ,
    int part,                 // Partition ID
    int K                     // Elements per partition
) {
    int partOffset = PART_OFFSET(part);
    
    // Temporary storage for reference derivatives
    real_t dUdr[K*NP], dUds[K*NP], dUdt[K*NP];
    
    // Apply differentiation matrices
    matMul_Dr_Partition(U, dUdr, part, K);
    matMul_Ds_Partition(U, dUds, part, K);
    matMul_Dt_Partition(U, dUdt, part, K);
    
    // Apply chain rule for physical derivatives
    for (int elem = 0; elem < K; ++elem) {
        for (int node = 0; node < NP; ++node) {
            int idx = partOffset + elem*NP + node;
            int geoIdx = idx * 9;  // 9 geometric factors per node
            
            DX[idx] = geoFactors[geoIdx+0]*dUdr[idx] + 
                      geoFactors[geoIdx+3]*dUds[idx] + 
                      geoFactors[geoIdx+6]*dUdt[idx];
                      
            DY[idx] = geoFactors[geoIdx+1]*dUdr[idx] + 
                      geoFactors[geoIdx+4]*dUds[idx] + 
                      geoFactors[geoIdx+7]*dUdt[idx];
                      
            DZ[idx] = geoFactors[geoIdx+2]*dUdr[idx] + 
                      geoFactors[geoIdx+5]*dUds[idx] + 
                      geoFactors[geoIdx+8]*dUdt[idx];
        }
    }
}
```

## Memory Management

### AllocateKernelMemory()

Allocates memory for all partitions in contiguous blocks:

```go
func (kp *KernelProgram) AllocateKernelMemory() error {
    floatSize := kp.GetFloatSize()
    intSize := kp.GetIntSize()
    
    // Total counts across all partitions
    totalNodes := kp.NumPartitions * kp.ElementsPerPart * kp.Np
    totalFaces := kp.NumPartitions * kp.ElementsPerPart * kp.Nfaces * kp.Nfp
    
    // Solution arrays - partition-blocked
    kp.memory["U"] = kp.device.Malloc(int64(totalNodes * floatSize), nil, nil)
    kp.memory["RHS"] = kp.device.Malloc(int64(totalNodes * floatSize), nil, nil)
    
    // Geometric factors - partition-blocked
    // Store as structure of arrays for better access patterns
    geoSize := totalNodes * 9 * floatSize  // 9 factors: rx,ry,rz,sx,sy,sz,tx,ty,tz
    kp.memory["geoFactors"] = kp.device.Malloc(int64(geoSize), nil, nil)
    
    // Face data - partition-blocked
    kp.memory["faceDataSend"] = kp.device.Malloc(int64(totalFaces * floatSize), nil, nil)
    kp.memory["faceDataRecv"] = kp.device.Malloc(int64(totalFaces * floatSize), nil, nil)
    
    // Partition connectivity for face exchange
    connSize := kp.NumPartitions * kp.ElementsPerPart * kp.Nfaces * 2  // (neighbor_part, neighbor_elem)
    kp.memory["partitionConn"] = kp.device.Malloc(int64(connSize * intSize), nil, nil)
    
    return nil
}
```

## Kernel Execution

### RunKernel Implementation

```go
func (kp *KernelProgram) RunKernel(name string, args ...interface{}) error {
    kernel, exists := kp.kernels[name]
    if !exists {
        return fmt.Errorf("kernel %s not found", name)
    }
    
    // Configure for partition-parallel execution
    // Each partition runs on its own GPU block / CPU thread
    outerDims := gocca.OCCADim{
        X: uint64(kp.NumPartitions),  // Partitions map to blocks
        Y: 1,
        Z: 1,
    }
    
    // Inner parallelism over nodes
    innerDims := gocca.OCCADim{
        X: uint64(kp.Np),  // Nodes map to threads
        Y: 1,
        Z: 1,
    }
    
    kernel.SetRunDims(outerDims, innerDims)
    
    return kernel.RunWithArgs(args...)
}
```

### Face Exchange Pattern

Inter-partition communication uses a separate kernel:

```c
@kernel void exchangeFaceData(
    const int Npart,
    const int K,
    const real_t* U,              // Solution data
    real_t* faceDataSend,         // Outgoing face data
    real_t* faceDataRecv,         // Incoming face data  
    const int_t* partitionConn,   // Connectivity info
    const int_t* faceNodeMap      // Maps volume nodes to face nodes
) {
    // Each partition prepares its boundary data
    for (int part = 0; part < Npart; ++part; @outer(0)) {
        
        // Extract face data for boundary elements
        for (int elem = 0; elem < K; ++elem) {
            for (int face = 0; face < NFACES; ++face) {
                
                // Check if this face is on partition boundary
                int connIdx = part*K*NFACES*2 + elem*NFACES*2 + face*2;
                int neighborPart = partitionConn[connIdx];
                
                if (neighborPart != part) {  // Inter-partition face
                    // Copy face data to send buffer
                    for (int fp = 0; fp < NFP; ++fp; @inner(0)) {
                        int volNode = faceNodeMap[face*NFP + fp];
                        int volIdx = ELEM_OFFSET(part, elem) + volNode;
                        int faceIdx = part*K*NFACES*NFP + elem*NFACES*NFP + face*NFP + fp;
                        
                        faceDataSend[faceIdx] = U[volIdx];
                    }
                }
            }
        }
    }
    
    // Actual communication happens between kernel calls
    // via CPU memcpy or GPU peer-to-peer transfers
}
```

## Example Kernels

### Complete RHS Computation

```c
@kernel void computeRHS(
    const int Npart,
    const int K,
    const real_t* U,
    const real_t* geoFactors,
    const real_t* faceDataRecv,
    real_t* RHS
) {
    // Partition-parallel execution
    for (int part = 0; part < Npart; ++part; @outer(0)) {
        
        // Allocate partition-local temporary arrays
        real_t DX[K*NP], DY[K*NP], DZ[K*NP];
        
        // Compute volume derivatives for this partition
        computePhysDerivPartition(U, geoFactors, DX, DY, DZ, part, K);
        
        // Process all elements in partition
        for (int elem = 0; elem < K; ++elem) {
            
            // Apply DG operations element-by-element
            for (int node = 0; node < NP; ++node; @inner(0)) {
                int idx = NODE_INDEX(part, elem, node);
                
                // Volume integral contribution
                real_t divF = computeDivergence(DX[idx], DY[idx], DZ[idx]);
                
                // Surface integral contribution
                real_t surfF = computeSurfaceFlux(elem, node, part, faceDataRecv);
                
                // Combine for RHS
                RHS[idx] = divF - surfF;
            }
        }
    }
}
```

## Performance Characteristics

### GPU Execution

For a typical problem on a modern GPU:
- **Global mesh**: 1,000,000 elements
- **GPU memory**: Can hold ~100,000 elements
- **Configuration**: Npart = 10, K = 100,000
- **Execution**:
   - 10 GPU blocks (one per partition)
   - Each block processes 100,000 elements
   - 35-100 threads per block (depending on polynomial order)
   - All partitions execute simultaneously

### Advantages

1. **Full GPU Utilization**: All streaming multiprocessors engaged
2. **Coalesced Memory Access**: Partition-blocked layout ensures adjacent threads access adjacent memory
3. **Minimal Synchronization**: Partitions execute independently
4. **Scalable**: Easy to adjust Npart based on problem size and hardware

### Design Decisions

1. **Why Partitions at @outer Level**:
   - Natural mapping to GPU blocks
   - Enables partition-independent execution
   - Simplifies face exchange pattern

2. **Why Element-Blocked Within Partitions**:
   - Maximizes cache reuse for element-local operations
   - Enables vectorization over nodes
   - Standard DG data layout

3. **Why Not Elements at @outer Level**:
   - Would require 1M blocks for 1M elements (exceeds GPU limits)
   - Partitioning provides right granularity for GPU blocks

## Usage Example

```go
// Configure for large 3D mesh
kp := NewKernelProgram(device, Config{
    NumPartitions:   64,      // 64 GPU blocks
    ElementsPerPart: 10000,   // 10K elements per partition
    Order:           4,       // 4th order polynomials
    FloatType:       Float64,
})

// Add DG operators
kp.AddStaticMatrix("Dr", Dr)
kp.AddStaticMatrix("Ds", Ds)
kp.AddStaticMatrix("Dt", Dt)

// Generate partition-aware code
preamble := kp.GenerateKernelMain()

// Build RHS kernel
kp.BuildKernel(rhsKernelSource, "computeRHS")

// Allocate memory for all partitions
kp.AllocateKernelMemory()

// Time stepping loop
for step := 0; step < nSteps; step++ {
    // Exchange face data between partitions
    kp.RunKernel("exchangeFaceData", ...)
    
    // Compute RHS with all partitions in parallel
    kp.RunKernel("computeRHS", 
        kp.NumPartitions,
        kp.ElementsPerPart,
        kp.GetMemory("U"),
        kp.GetMemory("geoFactors"),
        kp.GetMemory("faceDataRecv"),
        kp.GetMemory("RHS"))
    
    // Time integration
    kp.RunKernel("timeStep", ...)
}
```

## Conclusion

The partition-parallel execution model in KernelProgram provides:
- Efficient utilization of GPU/parallel hardware
- Natural mapping to OCCA's parallelism model
- Scalability from small to very large meshes
- Clear separation between independent partitions

This design enables high-performance DG computations while maintaining code clarity and portability across different parallel backends.
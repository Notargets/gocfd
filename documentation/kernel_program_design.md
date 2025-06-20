# KernelProgram Design Document

## Overview

The KernelProgram system provides a high-level abstraction for generating and executing GPU/accelerator kernels for Discontinuous Galerkin (DG) methods. It automates the generation of static data initialization, utility functions, and specialized operators while managing memory allocation and kernel execution through the gocca wrapper.

## Core Design Principles

1. **Element-Blocked Data Layout**: All arrays use element-blocked layout where nodes within an element are contiguous in memory. This is the universal DG pattern enabling optimal element-local operations.
2. **Partition-Based Parallelism**: Total elements divided into partitions, each executed as a single vectorized kernel
3. **Static Data Embedding**: Small matrices (Dr, Ds, Dt, LIFT) are compiled directly into kernel source as static arrays for optimal performance
4. **Vector Chaining**: Operations are designed to maximize throughput on long vectors (dimension K within each partition)
5. **Memory Persistence**: Device memory allocations persist across kernel executions
6. **Code Generation**: Kernel preambles containing static data and utilities are generated programmatically with layout-optimized operations

## Data Layout Convention

KernelProgram uses **element-blocked layout** throughout:
```
Array layout: [elem0_node0, elem0_node1, ..., elem0_nodeNP-1, elem1_node0, elem1_node1, ...]
```

This layout provides:
- **Cache efficiency**: All nodes of an element are accessed together during element-local operations
- **Vectorization**: Inner loops over nodes can be fully vectorized
- **Memory coalescing**: Adjacent threads access adjacent memory in GPU implementations
- **Natural DG structure**: Matches the mathematical formulation of DG methods

## Architecture

### Partitioning Model

The system divides the total number of elements (K_total) into partitions:
- Each partition contains K elements
- Each partition runs as a separate kernel execution or on a separate parallel unit
- Partitions communicate through face buffer exchanges
- Within each partition, operations are fully vectorized over element-blocked K×NP data

### KernelProgram Structure

```go
type KernelProgram struct {
    // Configuration
    Order           int                        // Polynomial order (N)
    Np              int                        // Number of nodes per element
    Nfp             int                        // Number of face points
    NumElements     int                        // Number of elements in this partition (K)
    
    // Data precision (default: Float64)
    FloatType       DataType                   // Float32 or Float64
    IntType         DataType                   // Int32 or Int64
    
    // Static data to embed
    StaticMatrices  map[string]utils.Matrix    // Dr, Ds, Dt, LIFT, etc.
    
    // Generated code
    kernelPreamble  string                     // Generated static data and utilities
    
    // Runtime resources
    device          *gocca.OCCADevice
    kernels         map[string]*gocca.OCCAKernel
    memory          map[string]*gocca.OCCAMemory
}
```

### Memory Layout

The system manages two types of data:

1. **Static Data** (compiled into kernels):
   - Differentiation matrices: Dr, Ds, Dt (Np × Np)
   - Lift matrix: LIFT (Np × 4*Nfp)
   - Small working arrays

2. **Dynamic Data** (allocated per partition with element-blocked layout):
   - Solution arrays: U[K][Np] - stored as contiguous array
   - Geometric factors: rx[K][Np], ry[K][Np], etc.
   - Face data buffers for inter-partition communication
   - Index arrays for connectivity

## Code Generation

### GenerateKernelMain()

This function generates the kernel preamble containing layout-optimized operations:

1. **Constants and Macros**
```c
#define ORDER 4
#define NP 35      // (N+1)*(N+2)*(N+3)/6
#define NFP 15     // (N+1)*(N+2)/2

// Type definitions based on precision setting
typedef double real_t;    // or float for Float32
typedef long   int_t;     // or int for Int32

// Element-blocked indexing helper
#define NODE_INDEX(elem, node) ((elem)*NP + (node))
```

2. **Layout-Optimized Matrix Operations**
```c
// Generated matrix multiplication for element-blocked data
// Optimized for cache efficiency with element-local operations
inline void matMul_Dr_Large(const real_t* U, real_t* result, int K) {
    for (int elem = 0; elem < K; ++elem) {
        for (int i = 0; i < NP; ++i) {
            real_t sum = REAL_ZERO;
            #pragma unroll
            for (int j = 0; j < NP; ++j) {
                sum += Dr[i][j] * U[elem*NP + j];
            }
            result[elem*NP + i] = sum;
        }
    }
}
```

The element-blocked layout enables:
- **Perfect cache reuse**: The static matrix Dr stays in cache while streaming through elements
- **Vectorization**: Modern compilers can vectorize the inner loop over nodes
- **Minimal memory movement**: Each element's data is accessed once in a contiguous block

3. **Specialized DG Operations**
```c
// Compute physical derivatives using element-blocked layout
inline void computePhysicalDerivatives(
    const real_t* U,         // Element-blocked: U[K][NP]
    const real_t* rx, const real_t* ry, const real_t* rz,
    const real_t* sx, const real_t* sy, const real_t* sz,
    const real_t* tx, const real_t* ty, const real_t* tz,
    real_t* DX, real_t* DY, real_t* DZ,
    int K
) {
    // Temporary arrays maintain element-blocked layout
    real_t dUdr[K*NP], dUds[K*NP], dUdt[K*NP];
    
    // Apply differentiation matrices
    matMul_Dr_Large(U, dUdr, K);
    matMul_Ds_Large(U, dUds, K);
    matMul_Dt_Large(U, dUdt, K);
    
    // Apply chain rule - perfect memory access pattern
    for (int elem = 0; elem < K; ++elem) {
        for (int node = 0; node < NP; ++node) {
            int idx = elem*NP + node;
            DX[idx] = rx[idx]*dUdr[idx] + sx[idx]*dUds[idx] + tx[idx]*dUdt[idx];
            DY[idx] = ry[idx]*dUdr[idx] + sy[idx]*dUds[idx] + ty[idx]*dUdt[idx];
            DZ[idx] = rz[idx]*dUdr[idx] + sz[idx]*dUds[idx] + tz[idx]*dUdt[idx];
        }
    }
}
```

## Memory Management

### AllocateKernelMemory()

Allocates persistent device memory for element-blocked runtime data:

```go
func (kp *KernelProgram) AllocateKernelMemory() error {
    floatSize := kp.GetFloatSize()
    intSize := kp.GetIntSize()
    
    // All arrays sized for element-blocked layout
    nodeCount := kp.Np * kp.NumElements  // Total nodes = Np × K
    
    // Solution arrays
    kp.memory["U"] = kp.device.Malloc(
        int64(nodeCount * floatSize),
        nil, 
        nil,
    )
    
    // Geometric factors maintain element-blocked layout
    geoFactors := []string{"rx", "ry", "rz", "sx", "sy", "sz", "tx", "ty", "tz"}
    for _, factor := range geoFactors {
        kp.memory[factor] = kp.device.Malloc(
            int64(nodeCount * floatSize),
            nil,
            nil,
        )
    }
    
    // Face data also element-blocked
    faceCount := 4 * kp.Nfp * kp.NumElements
    kp.memory["faceM"] = kp.device.Malloc(
        int64(faceCount * floatSize),
        nil,
        nil,
    )
    
    return nil
}
```

## Kernel Design Patterns

### Element-Local Operations

The element-blocked layout makes element-local operations extremely efficient:

```c
@kernel void elementLocalOperation(
    const int K,
    const real_t* U,      // U[K][NP]
    real_t* result        // result[K][NP]
) {
    for (int elem = 0; elem < K; ++elem; @outer) {
        // All nodes of element 'elem' are contiguous
        // Perfect for cache and vectorization
        for (int node = 0; node < NP; ++node; @inner) {
            int idx = elem*NP + node;
            result[idx] = someOperation(U[idx]);
        }
    }
}
```

### Face Operations

Face data extraction is straightforward with element-blocked layout:

```c
@kernel void extractFaceData(
    const int K,
    const real_t* U,           // U[K][NP]
    const int_t* faceNodes,    // Face node indices within element
    real_t* faceData          // faceData[K][Nfaces][Nfp]
) {
    for (int elem = 0; elem < K; ++elem; @outer) {
        for (int face = 0; face < 4; ++face) {
            for (int fp = 0; fp < NFP; ++fp; @inner) {
                int volumeNode = faceNodes[face*NFP + fp];
                int volumeIdx = elem*NP + volumeNode;
                int faceIdx = elem*4*NFP + face*NFP + fp;
                faceData[faceIdx] = U[volumeIdx];
            }
        }
    }
}
```

## Performance Benefits

The element-blocked layout provides significant performance advantages:

1. **Cache Efficiency**:
   - Static matrices remain in L1/L2 cache
   - Element data accessed in contiguous blocks
   - Minimal cache misses during element-local operations

2. **Vectorization**:
   - Inner loops over NP nodes can be fully vectorized
   - Compiler can optimize for SIMD instructions
   - Natural alignment for vector loads/stores

3. **Memory Bandwidth**:
   - Contiguous memory access maximizes bandwidth utilization
   - Coalesced memory access on GPUs
   - Reduced memory latency

4. **Parallelization**:
   - Elements are independent computational units
   - Perfect for GPU thread blocks or CPU SIMD lanes
   - No false sharing between elements

## Usage Example

```go
func main() {
    // Create KernelProgram with element-blocked layout
    device, _ := gocca.NewDevice(`{"mode": "CUDA"}`)
    kp := NewKernelProgram(device, Config{
        Order:       4,
        NumElements: 1000,  // K elements in this partition
        FloatType:   Float64,
    })
    
    // Add DG operators
    kp.AddStaticMatrix("Dr", Dr)
    kp.AddStaticMatrix("Ds", Ds)
    kp.AddStaticMatrix("Dt", Dt)
    
    // Generate optimized code for element-blocked layout
    preamble := kp.GenerateKernelMain()
    
    // User kernels work with element-blocked data
    rhsKernel := `
    @kernel void computeRHS(
        const int K,
        const real_t* U,      // Element-blocked
        real_t* RHS           // Element-blocked
    ) {
        // Temporary storage
        real_t DX[NP*K], DY[NP*K], DZ[NP*K];
        
        // Compute derivatives - optimized for element-blocked layout
        computePhysicalDerivatives(U, rx, ry, rz, sx, sy, sz, tx, ty, tz,
                                  DX, DY, DZ, K);
        
        // Further operations maintain element-blocked pattern
        for (int elem = 0; elem < K; ++elem) {
            for (int node = 0; node < NP; ++node) {
                int idx = elem*NP + node;
                RHS[idx] = /* flux divergence using DX, DY, DZ */;
            }
        }
    }
    `
    
    // Build and execute
    kp.BuildKernel(rhsKernel, "computeRHS")
    kp.AllocateKernelMemory()
    
    // Initialize with element-blocked data
    hostU := make([]float64, kp.Np * kp.NumElements)
    // Fill hostU with element-blocked layout
    kp.GetMemory("U").CopyFrom(unsafe.Pointer(&hostU[0]))
    
    // Execute optimized kernel
    kp.RunKernel("computeRHS", kp.NumElements, 
                 kp.GetMemory("U"), kp.GetMemory("RHS"))
}
```

## Conclusion

The KernelProgram system's commitment to element-blocked data layout enables the generation of highly optimized kernels for DG methods. By standardizing on this layout pattern, the system can:

- Generate cache-optimal matrix operations
- Enable aggressive compiler optimizations
- Maximize hardware utilization
- Simplify kernel development

This design choice reflects the reality that all production DG codes use element-blocked layout, making it the natural choice for a DG-focused kernel generation system.
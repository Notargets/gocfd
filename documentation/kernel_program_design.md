# KernelProgram Design Document

## Overview

The KernelProgram system provides a high-level abstraction for generating and executing GPU/accelerator kernels for Discontinuous Galerkin (DG) methods. It automates the generation of static data initialization, utility functions, and specialized operators while managing memory allocation and kernel execution through the gocca wrapper.

## Core Design Principles

1. **Partition-Based Parallelism**: Total elements divided into partitions, each executed as a single vectorized kernel
2. **Static Data Embedding**: Small matrices (Dr, Ds, Dt, LIFT) are compiled directly into kernel source as static arrays for optimal performance
3. **Vector Chaining**: Operations are designed to maximize throughput on long vectors (dimension K within each partition)
4. **Memory Persistence**: Device memory allocations persist across kernel executions
5. **Code Generation**: Kernel preambles containing static data and utilities are generated programmatically

## Architecture

### Partitioning Model

The system divides the total number of elements (K_total) into partitions:
- Each partition contains K elements
- Each partition runs as a separate kernel execution or on a separate parallel unit
- Partitions communicate through face buffer exchanges
- Within each partition, operations are fully vectorized over K×NP dimensions

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

type DataType int
const (
    Float32 DataType = iota
    Float64
    Int32
    Int64
)
```

### Memory Layout

The system manages two types of data:

1. **Static Data** (compiled into kernels):
    - Differentiation matrices: Dr, Ds, Dt (Np × Np)
    - Lift matrix: LIFT (Np × 4*Nfp)
    - Small working arrays

2. **Dynamic Data** (allocated per partition):
    - Solution arrays: U (Np × K) - entire partition
    - Geometric factors: rx, ry, rz, sx, sy, sz, tx, ty, tz (Np × K) - entire partition
    - Face data buffers for inter-partition communication
    - Index arrays for connectivity

Key principle: Operations treat (Np × K) arrays as single vectorized entities, maximizing throughput on vector architectures.

## Code Generation

### GenerateKernelMain()

This function generates the kernel preamble containing:

1. **Constants and Macros**
```c
#define ORDER 4
#define NP 35      // (N+1)*(N+2)*(N+3)/6
#define NFP 15     // (N+1)*(N+2)/2

// Type definitions based on precision setting
typedef double real_t;    // or float for Float32
typedef long   int_t;     // or int for Int32
```

2. **Static Matrix Initialization** (type based on FloatType)
```c
__constant__ real_t Dr[NP][NP] = {  // real_t defined as float or double
    {1.234, 5.678, ...},
    ...
};
```

3. **Utility Functions**
```c
// Matrix-vector multiplication using static Dr
inline void applyDr(const real_t* u, real_t* du) {
    for (int i = 0; i < NP; ++i) {
        real_t sum = 0.0;
        #pragma unroll
        for (int j = 0; j < NP; ++j) {
            sum += Dr[i][j] * u[j];
        }
        du[i] = sum;
    }
}
```

4. **Specialized Operators**
```c
// Fused derivative chain for entire partition
inline void computePhysicalDerivativesPartition(
    const real_t* U,         // (NP, K) input
    const real_t* rx,        // (NP, K) geometric factors  
    const real_t* ry,
    const real_t* rz,
    const real_t* sx,
    const real_t* sy,
    const real_t* sz,
    const real_t* tx,
    const real_t* ty,
    const real_t* tz,
    real_t* DX,              // (NP, K) output
    real_t* DY,
    real_t* DZ,
    int K
) {
    // Temporary arrays for reference derivatives
    real_t dUdr[NP*K];
    real_t dUds[NP*K];
    real_t dUdt[NP*K];
    
    // Apply differentiation matrices
    matMul_Dr_Large(U, dUdr, K);
    matMul_Ds_Large(U, dUds, K);
    matMul_Dt_Large(U, dUdt, K);
    
    // Apply chain rule - fully vectorized
    for (int idx = 0; idx < NP*K; ++idx) {
        DX[idx] = rx[idx]*dUdr[idx] + sx[idx]*dUds[idx] + tx[idx]*dUdt[idx];
        DY[idx] = ry[idx]*dUdr[idx] + sy[idx]*dUds[idx] + ty[idx]*dUdt[idx];
        DZ[idx] = rz[idx]*dUdr[idx] + sz[idx]*dUds[idx] + tz[idx]*dUdt[idx];
    }
}
```

## Memory Management

### AllocateKernelMemory()

Allocates persistent device memory for runtime data:

```go
func (kp *KernelProgram) AllocateKernelMemory() error {
    // Get size based on precision
    floatSize := 8  // Float64
    if kp.FloatType == Float32 {
        floatSize = 4
    }
    
    // Solution arrays
    nodeCount := kp.Np * kp.NumElements
    kp.memory["U"] = kp.device.Malloc(
        int64(nodeCount * floatSize),
        nil, 
        nil,
    )
    
    // Geometric factors (stored separately for vectorization)
    geoFactors := []string{"rx", "ry", "rz", "sx", "sy", "sz", "tx", "ty", "tz"}
    for _, factor := range geoFactors {
        kp.memory[factor] = kp.device.Malloc(
            int64(nodeCount * floatSize),
            nil,
            nil,
        )
    }
    
    // Face buffers, index arrays, etc.
    // ...
    
    return nil
}
```

## Kernel Execution

### Monolithic Kernel Design

The system favors large, monolithic kernels that perform complete computational phases without intermediate memory transfers. This approach:

1. **Minimizes kernel launch overhead** - One launch instead of many
2. **Maximizes data reuse** - Intermediate results stay in registers/fast memory
3. **Reduces memory bandwidth** - No need to write/read intermediate arrays
4. **Enables cross-operation optimization** - Compiler can optimize across operation boundaries

Occa supports arbitrarily complex kernels (within hardware instruction limits, typically millions of instructions), making it ideal for combining multiple DG operations into single kernels.

### Generated Kernels

The system supports generating various computational kernels:

1. **Volume Derivatives Kernel** (operates on entire partition)
```c
@kernel void computeVolumeDerivatives(
    const int K,              // Number of elements in this partition
    const real_t* U,          @ restrict  // (NP, K) - entire partition
    const real_t* rx,         @ restrict  // (NP, K) geometric factors
    const real_t* ry,         @ restrict
    const real_t* rz,         @ restrict
    const real_t* sx,         @ restrict
    const real_t* sy,         @ restrict
    const real_t* sz,         @ restrict
    const real_t* tx,         @ restrict
    const real_t* ty,         @ restrict
    const real_t* tz,         @ restrict
    real_t* DX,               @ restrict  // (NP, K) output
    real_t* DY,               @ restrict
    real_t* DZ                @ restrict
) {
    // Temporary arrays for reference derivatives - entire partition
    real_t dUdr[NP*K];
    real_t dUds[NP*K];
    real_t dUdt[NP*K];
    
    // Apply differentiation matrices to entire partition
    // Dr @ U where Dr is (NP×NP) and U is (NP×K)
    matMul_Dr_Large(U, dUdr, K);  // Vectorized over entire K
    matMul_Ds_Large(U, dUds, K);
    matMul_Dt_Large(U, dUdt, K);
    
    // Apply chain rule - vectorized over entire NP*K array
    for (int idx = 0; idx < NP*K; ++idx) {
        DX[idx] = rx[idx]*dUdr[idx] + sx[idx]*dUds[idx] + tx[idx]*dUdt[idx];
        DY[idx] = ry[idx]*dUdr[idx] + sy[idx]*dUds[idx] + ty[idx]*dUdt[idx];
        DZ[idx] = rz[idx]*dUdr[idx] + sz[idx]*dUds[idx] + tz[idx]*dUdt[idx];
    }
}

// Optimized matrix multiply for large K
inline void matMul_Dr_Large(const real_t* U, real_t* result, int K) {
    // This entire operation is vectorized
    // Modern compilers/accelerators can vectorize across K dimension
    for (int i = 0; i < NP; ++i) {
        for (int k = 0; k < K; ++k) {
            real_t sum = 0.0;
            #pragma unroll
            for (int j = 0; j < NP; ++j) {
                sum += Dr[i][j] * U[j*K + k];
            }
            result[i*K + k] = sum;
        }
    }
}
```

2. **Surface Integral Kernel**
```c
@kernel void applySurfaceIntegral(
    const int K,
    const float* faceFlux,   @ restrict  // (4*NFP, K)
    float* rhsU              @ restrict  // (NP, K)
) {
    for (int elem = 0; elem < K; ++elem; @outer(0)) {
        @shared float localFlux[4*NFP];
        
        // Load face fluxes
        for (int i = 0; i < 4*NFP; ++i; @inner(0)) {
            if (i < 4*NFP) {
                localFlux[i] = faceFlux[i + elem*4*NFP];
            }
        }
        
        @barrier("local");
        
        // Apply LIFT matrix
        for (int i = 0; i < NP; ++i; @inner(0)) {
            if (i < NP) {
                float sum = 0.0f;
                #pragma unroll
                for (int j = 0; j < 4*NFP; ++j) {
                    sum += LIFT[i][j] * localFlux[j];
                }
                rhsU[i + elem*NP] += sum;
            }
        }
    }
}
```

## Example Usage

```go
package main

import (
    "github.com/notargets/gocca"
    "myproject/kernelprogram"
    "myproject/utils"
)

func main() {
    // 1. Initialize device
    device, err := gocca.NewDevice(`{"mode": "CUDA", "device_id": 0}`)
    if err != nil {
        panic(err)
    }
    defer device.Free()
    
    // 2. Create KernelProgram
    kp := kernelprogram.New(device, kernelprogram.Config{
        Order:       4,
        NumElements: 10000,      // Elements in this partition
        FloatType:   kernelprogram.Float64,  // Optional, defaults to Float64
        IntType:     kernelprogram.Int64,    // Optional, defaults to Int64
    })
    
    // 3. Add static matrices
    Dr, Ds, Dt := computeDerivativeMatrices(kp.Order)
    LIFT := computeLiftMatrix(kp.Order)
    
    kp.AddStaticMatrix("Dr", Dr)
    kp.AddStaticMatrix("Ds", Ds)
    kp.AddStaticMatrix("Dt", Dt)
    kp.AddStaticMatrix("LIFT", LIFT)
    
    // 4. Generate kernel preamble
    preamble := kp.GenerateKernelMain()
    
    // 5. Build kernels with preamble + handwritten source
    volumeKernelSource := preamble + volumeKernelCode
    volumeKernel, err := device.BuildKernelFromString(
        volumeKernelSource, 
        "computeVolumeDerivatives",
        nil,
    )
    if err != nil {
        panic(err)
    }
    kp.RegisterKernel("volumeDerivatives", volumeKernel)
    
    // 6. Allocate memory
    err = kp.AllocateKernelMemory()
    if err != nil {
        panic(err)
    }
    
    // 7. Initialize solution data
    U := kp.GetMemory("U")
    U.CopyFrom(initialSolution)
    
    // 8. Execute monolithic RHS kernel
    err = kp.RunKernel("computeRHS", 
        kp.NumElements,  // K for this partition
        kp.GetMemory("U"),
        kp.GetMemory("faceM"),
        kp.GetMemory("faceP"),
        kp.GetMemory("faceTypes"),
        kp.GetMemory("bcData"),
        kp.GetMemory("rx"),
        kp.GetMemory("ry"),
        kp.GetMemory("rz"),
        kp.GetMemory("sx"),
        kp.GetMemory("sy"),
        kp.GetMemory("sz"),
        kp.GetMemory("tx"),
        kp.GetMemory("ty"),
        kp.GetMemory("tz"),
        kp.GetMemory("nx"),
        kp.GetMemory("ny"),
        kp.GetMemory("nz"),
        kp.GetMemory("Fscale"),
        kp.GetMemory("RHS"),
    )
    if err != nil {
        panic(err)
    }
    
    // 9. Retrieve results
    hostDX := make([]float64, kp.Np * kp.NumElements)
    kp.GetMemory("DX").CopyTo(unsafe.Pointer(&hostDX[0]))
}
```

## Implementation Details

### Static Matrix Formatting

The `formatStaticMatrix` function converts utils.Matrix to C array initialization:

```go
func (kp *KernelProgram) formatStaticMatrix(name string, m utils.Matrix) string {
    rows, cols := m.Dims()
    var sb strings.Builder
    
    // Use appropriate type based on FloatType
    typeStr := "double"
    if kp.FloatType == Float32 {
        typeStr = "float"
    }
    
    sb.WriteString(fmt.Sprintf("__constant__ %s %s[%d][%d] = {\n", 
        typeStr, name, rows, cols))
    
    for i := 0; i < rows; i++ {
        sb.WriteString("    {")
        for j := 0; j < cols; j++ {
            if j > 0 {
                sb.WriteString(", ")
            }
            // Use appropriate precision suffix
            if kp.FloatType == Float32 {
                sb.WriteString(fmt.Sprintf("%.7ef", m.At(i, j)))
            } else {
                sb.WriteString(fmt.Sprintf("%.15e", m.At(i, j)))
            }
        }
        sb.WriteString("},\n")
    }
    sb.WriteString("};\n\n")
    
    return sb.String()
}
```

### Kernel Execution Wrapper

```go
func (kp *KernelProgram) RunKernel(name string, args ...interface{}) error {
    kernel, exists := kp.kernels[name]
    if !exists {
        return fmt.Errorf("kernel %s not found", name)
    }
    
    // Single kernel execution for entire partition
    // No thread dimensions needed - kernel internally vectorizes
    kernel.SetRunDims(
        gocca.OCCADim{X: 1, Y: 1, Z: 1},  // Single execution
        gocca.OCCADim{X: 1, Y: 1, Z: 1},
    )
    
    return kernel.RunWithArgs(args...)
}
```

## Performance Considerations

1. **Static Data**: Compiled into kernel, likely placed in constant cache
2. **Vectorization**: Operations on K elements provide maximum vector length for hardware units
3. **Memory Bandwidth**: Streaming entire (NP×K) arrays maximizes bandwidth utilization
4. **Cache Reuse**: Small static matrices (NP×NP) stay in cache while long vectors stream through
5. **Kernel Fusion**: Monolithic kernels keep intermediate results in fast memory, reducing bandwidth requirements
6. **Partition Communication**: Face buffers enable overlap of computation and communication
7. **Precision Selection**:
    - Float32: Doubles bandwidth and compute throughput on most hardware
    - Float64: Required for stability in some DG applications
    - Mixed precision possible: Float32 for bandwidth-limited ops, Float64 for accuracy-critical sections

## Future Extensions

1. **Multi-GPU Support**: Each partition on a separate GPU with face buffer exchanges via MPI/NCCL
2. **Mixed Precision**: Use float32 for bandwidth-limited operations
3. **Kernel Caching**: Save compiled binaries for faster startup
4. **Auto-tuning**: Profile and select optimal partition sizes
5. **Template Specialization**: Generate kernels for specific polynomial orders

## Multi-Partition Execution Model

In a full implementation with multiple partitions:

```go
// Pseudo-code for multi-partition execution
type PartitionedSolver struct {
    partitions []KernelProgram
    faceBuffers []FaceBuffer
}

func (ps *PartitionedSolver) TimeStep() {
    // 1. Exchange face data between partitions
    // (Could be overlapped with previous timestep computation)
    ps.ExchangeFaceData()
    
    // 2. Each partition computes complete RHS in one kernel
    for _, partition := range ps.partitions {
        partition.RunKernel("computeRHS", 
            // All required memory arguments
        )
    }
    
    // 3. Time integration (could also be fused into RHS kernel)
    for _, partition := range ps.partitions {
        partition.RunKernel("timeIntegration", 
            partition.GetMemory("U"),
            partition.GetMemory("RHS"),
            dt,
        )
    }
}
```

This design minimizes synchronization points and maximizes the work done in each kernel launch.
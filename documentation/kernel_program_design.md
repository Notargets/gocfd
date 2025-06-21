# KernelProgram Design Document

## Overview

The KernelProgram system provides infrastructure for managing GPU/accelerator kernels with partition-parallel execution. It handles memory allocation, kernel compilation, and code generation to simplify writing high-performance parallel kernels. The system is designed as a general-purpose tool that users configure for their specific computational needs.

## Core Design Principles

1. **User-Controlled Memory Layout**: Users specify what arrays to allocate, their sizes, and alignment requirements
2. **Partition-Parallel Support**: Handles variable partition sizes through K[part] array
3. **Automatic Offset Management**: Generates offset arrays for easy partition-based indexing in kernels
4. **Static Matrix Embedding**: Embeds matrices as compile-time constants for optimal performance
5. **Vectorizable Matrix Operations**: Provides macros for matrix multiplication that compilers can vectorize
6. **OCCA Compliance**: Works with all OCCA backends maintaining proper @outer/@inner structure

## Core Structure

```go
type KernelProgram struct {
    // Partition configuration
    NumPartitions int
    K []int  // Variable elements per partition
    
    // Type configuration
    FloatType, IntType DataType
    
    // Static data for embedding
    StaticMatrices map[string]Matrix
    
    // Array tracking for macro generation
    allocatedArrays []string
    
    // Generated code
    kernelPreamble string
    
    // Runtime resources
    device *gocca.OCCADevice
    kernels map[string]*gocca.OCCAKernel
    pooledMemory map[string]*gocca.OCCAMemory
}
```

## Memory Management

### User-Specified Array Allocation

Users define their memory requirements through ArraySpec:

```go
type ArraySpec struct {
    Name      string
    Size      int64          // Total size in bytes
    Alignment AlignmentType  // Alignment requirement
}

type AlignmentType int
const (
    NoAlignment      AlignmentType = 1     
    CacheLineAlign   AlignmentType = 64    
    WarpAlign        AlignmentType = 128   
    PageAlign        AlignmentType = 4096  
)
```

### Automatic Offset Generation

For each allocated array, KernelProgram automatically generates offset arrays and access macros:

```go
func (kp *KernelProgram) AllocateArrays(specs []ArraySpec) error {
    for _, spec := range specs {
        // Allocate aligned memory
        memory := device.Malloc(spec.Size, spec.Alignment)
        kp.pooledMemory[spec.Name + "_global"] = memory
        
        // Generate per-partition offsets
        offsets := kp.calculateOffsets(spec)
        offsetMem := device.Malloc(len(offsets) * 8)
        offsetMem.CopyFrom(offsets)
        kp.pooledMemory[spec.Name + "_offsets"] = offsetMem
        
        // Track array name for macro generation
        kp.allocatedArrays = append(kp.allocatedArrays, spec.Name)
    }
}
```

### Kernel Access Pattern

Kernels access partition data using generated macros:

```c
@kernel void myKernel(
    const int* K,
    const real_t* U_global,
    const int_t* U_offsets,
    real_t* RHS_global,
    const int_t* RHS_offsets
) {
    for (int part = 0; part < NPART; ++part; @outer) {
        for (int node = 0; node < NP; ++node; @inner) {
            // Access partition data through macros
            const real_t* U = U_g(part);
            real_t* RHS = RHS_g(part);
            int K_part = K[part];
            
            // Process elements in this partition
            for (int elem = 0; elem < K_part; ++elem) {
                // Work with U[elem * NP + node]
                RHS[elem * NP + node] = 2.0 * U[elem * NP + node];
            }
        }
    }
}
```

## Code Generation

### Kernel Preamble

The generated preamble provides essential definitions and utilities:

```c
// Type definitions
typedef double real_t;  // or float based on configuration
typedef int64_t int_t;  // or int32_t based on configuration

// Constants (example values)
#define NP 20       // Nodes per element
#define NFP 10      // Face nodes per element  
#define NFACES 4    // Faces per element
#define NPART 64    // Number of partitions

// Static matrices embedded as constants
const real_t Dr[20][20] = { /* values */ };
const real_t Ds[20][20] = { /* values */ };
const real_t LIFT[20][40] = { /* values */ };

// Vectorizable matrix multiplication macros
#define MATMUL_Dr(IN, OUT, K) do { \
    for (int i = 0; i < 20; ++i) { \
        for (int elem = 0; elem < (K); ++elem) { \
            real_t sum = 0.0; \
            for (int j = 0; j < 20; ++j) { \
                sum += Dr[i][j] * (IN)[elem * 20 + j]; \
            } \
            (OUT)[elem * 20 + i] = sum; \
        } \
    } \
} while(0)

// Similar macros generated for each static matrix

// Partition access macros (generated for each allocated array)
#define U_PART(part) (U_global + U_offsets[part])
#define RHS_PART(part) (RHS_global + RHS_offsets[part])
#define rx_PART(part) (rx_global + rx_offsets[part])
// ... similar macros for each user-specified array
```

### Partition Access Macros

For each allocated array, KernelProgram generates a macro for easy partition access:

```c
// For array "U", generates:
#define U_PART(part) (U_global + U_offsets[part])

// Usage in kernel:
const real_t* U = U_PART(part);
```

This approach:
- Eliminates manual offset calculations
- Makes kernels cleaner and more readable
- Maintains zero overhead (macros expand at compile time)
- Prevents indexing errors

### Naming Convention

For each array specified by the user, KernelProgram generates:
- Global array parameter: `<name>_global`
- Offset array parameter: `<name>_offsets`
- Access macro: `<name>_PART(part)`

Example: Array "U" becomes `U_global`, `U_offsets`, and `U_PART(part)`

### Matrix Operation Macros

Two categories of matrix multiplication macros are generated:

1. **Square Matrix Operations** (e.g., differentiation matrices):
```c
#define MATMUL_Dr(IN, OUT, K) do { \
    for (int i = 0; i < NP; ++i) { \
        for (int elem = 0; elem < (K); ++elem) { \
            real_t sum = 0.0; \
            for (int j = 0; j < NP; ++j) { \
                sum += Dr[i][j] * (IN)[elem * NP + j]; \
            } \
            (OUT)[elem * NP + i] = sum; \
        } \
    } \
} while(0)
```

2. **Rectangular Matrix Operations** (e.g., LIFT matrix):
```c
#define MATMUL_LIFT(IN, OUT, K) do { \
    for (int i = 0; i < NP; ++i) { \
        for (int elem = 0; elem < (K); ++elem) { \
            real_t sum = 0.0; \
            for (int j = 0; j < NFACES * NFP; ++j) { \
                sum += LIFT[i][j] * (IN)[elem * NFACES * NFP + j]; \
            } \
            (OUT)[elem * NP + i] = sum; \
        } \
    } \
} while(0)
```

These macros:
- Expose all loops to the compiler for vectorization
- Work with variable K values per partition
- Avoid function call overhead
- Enable architecture-specific optimizations

## Kernel Execution

KernelProgram manages parameter passing with the renamed arrays:

```go
func (kp *KernelProgram) RunKernel(name string, args ...interface{}) error {
    kernel, exists := kp.kernels[name]
    if !exists {
        return fmt.Errorf("kernel %s not found", name)
    }
    
    // KernelProgram automatically handles the _global/_offsets naming
    // when building the argument list
    return kernel.RunWithArgs(args...)
}

// Example: When user calls with "U", "RHS" arrays, KernelProgram passes:
// U_global, U_offsets, RHS_global, RHS_offsets to the kernel
```

## Usage Example

### Application Setup

```go
// Configure memory layout for a DG solver
floatSize := 8  // sizeof(double)
intSize := 8    // sizeof(int64)

// Define memory requirements
arrays := []ArraySpec{
    // Solution arrays
    {Name: "U",    Size: totalNodes * floatSize, Alignment: CacheLineAlign},
    {Name: "RHS",  Size: totalNodes * floatSize, Alignment: CacheLineAlign},
    
    // Geometric factors  
    {Name: "rx",   Size: totalNodes * floatSize, Alignment: NoAlignment},
    {Name: "ry",   Size: totalNodes * floatSize, Alignment: NoAlignment},
    
    // Face data
    {Name: "nx",   Size: totalFaces * floatSize, Alignment: NoAlignment},
    {Name: "ny",   Size: totalFaces * floatSize, Alignment: NoAlignment},
}

// Create kernel program
kp := NewKernelProgram(device, Config{
    NumPartitions: nPart,
    K:            partitionSizes,  // From mesh partitioner
    FloatType:    Float64,
    IntType:      Int64,
})

// Add static matrices
kp.AddStaticMatrix("Dr", Dr)
kp.AddStaticMatrix("Ds", Ds)
kp.AddStaticMatrix("LIFT", LIFT)

// Generate preamble with macros
kp.GeneratePreamble()

// Allocate memory with automatic offset generation
kp.AllocateArrays(arrays)

// Build kernels
kp.BuildKernel(volumeKernelSource, "volumeKernel")
kp.BuildKernel(surfaceKernelSource, "surfaceKernel")
```

### Kernel Implementation

```c
@kernel void volumeKernel(
    const int* K,
    const real_t* U_global,     const int_t* U_offsets,
    const real_t* rx_global,    const int_t* rx_offsets,
    const real_t* ry_global,    const int_t* ry_offsets,
    real_t* RHS_global,         const int_t* RHS_offsets
) {
    for (int part = 0; part < NPART; ++part; @outer) {
        for (int node = 0; node < NP; ++node; @inner) {
            const int K_part = K[part];
            
            // Get partition pointers using macros
            const real_t* U = U_g(part);
            const real_t* rx = rx_g(part);
            const real_t* ry = ry_g(part);
            real_t* RHS = RHS_g(part);
            
            // Temporary arrays for derivatives
            real_t Ur[K_part * NP];
            real_t Us[K_part * NP];
            
            // Use vectorizable macros
            MATMUL_Dr(U, Ur, K_part);
            MATMUL_Ds(U, Us, K_part);
            
            // Apply chain rule
            for (int elem = 0; elem < K_part; ++elem) {
                int idx = elem * NP + node;
                RHS[idx] = rx[idx] * Ur[idx] + ry[idx] * Us[idx];
            }
        }
    }
}
```

## Design Benefits

1. **Flexibility**: Users control memory layout for their specific applications
2. **Performance**: Vectorizable macros enable compiler optimizations
3. **Simplicity**: Automatic offset generation and partition access macros simplify kernel code
4. **Portability**: Works across all OCCA backends without modification
5. **Efficiency**: Static matrix embedding eliminates memory accesses
6. **Scalability**: Variable partition sizes support real-world mesh partitioning
7. **Safety**: Partition access macros prevent manual indexing errors
8. **Clarity**: Clean kernel code with minimal boilerplate

This design provides a clean infrastructure layer that handles the complexity of memory management and code generation while giving users full control over their computational patterns.
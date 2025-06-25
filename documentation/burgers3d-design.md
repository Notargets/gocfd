# 3D Inviscid Burger’s Equation Solver Implementation Guide

## Overview

This document describes the implementation of a 3D scalar inviscid Burger’s equation solver using discontinuous Galerkin (DG) methods in conservation form. The solver combines `kernel_program.go` from gocca with the tetrahedra package from gocfd to create a partition-parallel solver supporting both OpenMP and CUDA execution.

The scalar Burger’s equation serves as an ideal test problem for the numerical methods and infrastructure that will be used in subsequent Navier-Stokes solvers, providing nonlinear flux terms and shock formation while maintaining simplicity.

## Mathematical Formulation

The 3D scalar inviscid Burger’s equation in conservation form:

```
∂u/∂t + ∂F/∂x + ∂G/∂y + ∂H/∂z = 0

where:
F = u²/2
G = u²/2  
H = u²/2
```

In the DG formulation:

```
∂u/∂t = -∇·F^v + LIFT * (F* - F)·n̂

where:
- F^v is the volume flux contribution
- F* is the numerical flux at faces
- F is the local flux at faces
- n̂ is the outward normal
```

## Architecture Components

### 1. Driver Program Structure

The driver program (`burgers3d.go`) resides in the gocfd repository and coordinates:

- Mesh reading and preprocessing
- Element matrix setup via DG3D
- Device initialization (OpenMP/CUDA)
- Kernel program creation and execution
- Time stepping loop

### 2. Data Flow

```
Mesh File → Element3D → SplitMesh → KernelProgram → Time Integration
             ↓                         ↓
         Geometric Factors         Static Matrices
                                   Face Buffers
```

## Implementation Steps

### Step 1: Mesh Reading and Element Setup

```go
// Read mesh using DG3D mesh reader
el, err := tetelement.NewElement3D(order, meshFile)

// This single call:
// - Reads the mesh file
// - Computes reference element operators (Dr, Ds, Dt)
// - Calculates physical coordinates (X, Y, Z)
// - Computes geometric factors (Rx,Ry,Rz,Sx,Sy,Sz,Tx,Ty,Tz,J)
// - Builds connectivity maps (VmapM, VmapP)
// - Sets up boundary condition mappings
// - Splits mesh into partitions if EToP is present
```

### Step 2: Device Setup

```go
// Initialize gocca device
var device *gocca.OCCADevice
if useCUDA {
    device = gocca.DeviceFromString("mode: 'CUDA', device_id: 0")
} else {
    device = gocca.DeviceFromString("mode: 'OpenMP'")
}
defer device.Free()
```

### Step 3: Extract Partition Information

```go
// Handle both partitioned and non-partitioned meshes
var numPartitions int
var K []int
var workingElements []*tetelement.Element3D

if el.Split != nil {
    // Partitioned mesh
    numPartitions = len(el.Split)
    K = make([]int, numPartitions)
    workingElements = el.Split
    for i, partEl := range el.Split {
        K[i] = partEl.K
    }
} else {
    // Non-partitioned mesh - treat as single partition
    numPartitions = 1
    K = []int{el.K}
    workingElements = []*tetelement.Element3D{el}
}

// Verify CUDA limits
if device.Mode() == "CUDA" {
    kpartMax := 0
    for _, k := range K {
        if k > kpartMax {
            kpartMax = k
        }
    }
    if kpartMax > 1024 {
        log.Fatal("CUDA @inner limit exceeded: max partition size is 1024")
    }
}
```

### Step 4: Create KernelProgram with Static Matrices

```go
// Initialize kernel program
kp := kernel_program.NewKernelProgram(device, kernel_program.Config{
    K:         K,
    FloatType: kernel_program.Float64,
    IntType:   kernel_program.INT64,
})

// Add static matrices from Element3D (reference element)
// These are small Np x Np matrices shared across all elements
kp.AddStaticMatrix("Dr", el.Dr)
kp.AddStaticMatrix("Ds", el.Ds)
kp.AddStaticMatrix("Dt", el.Dt)
kp.AddStaticMatrix("MassMatrix", el.MassMatrix)
kp.AddStaticMatrix("LIFT", el.LIFT)
```

### Step 5: Calculate Memory Requirements and Allocate Arrays

The derivative arrays (ur, us, ut, etc.) and flux arrays are allocated as device memory to avoid stack overflow issues with large temporary arrays and enable reuse across multiple kernel calls.

```go
// Calculate sizes
Np := el.Np
Nfp := el.Nfp
totalSolutionNodes := 0
totalFaceNodes := 0
totalGeometricFactors := 0

for _, k := range K {
    totalSolutionNodes += Np * k
    totalFaceNodes += Nfp * 4 * k  // 4 faces per tet
    totalGeometricFactors += Np * k
}

// Define array specifications
arrays := []kernel_program.ArraySpec{
    // Solution array
    {
        Name:      "u",
        Size:      int64(totalSolutionNodes * 8), // 8 bytes per float64
        DataType:  kernel_program.Float64,
        Alignment: kernel_program.CacheLineAlign,
    },
    // M buffer for face data
    {
        Name:      "M",
        Size:      int64(totalFaceNodes * 8),
        DataType:  kernel_program.Float64,
        Alignment: kernel_program.CacheLineAlign,
    },
    // Geometric factors (10 total: Rx,Ry,Rz,Sx,Sy,Sz,Tx,Ty,Tz,J)
    {
        Name:      "Rx",
        Size:      int64(totalGeometricFactors * 8),
        DataType:  kernel_program.Float64,
        Alignment: kernel_program.CacheLineAlign,
    },
    {
        Name:      "Ry",
        Size:      int64(totalGeometricFactors * 8),
        DataType:  kernel_program.Float64,
        Alignment: kernel_program.CacheLineAlign,
    },
    // ... repeat for Rz, Sx, Sy, Sz, Tx, Ty, Tz, J
    
    // RHS arrays for RK5SSP4 (5 stages)
    {
        Name:      "rhs0",
        Size:      int64(totalSolutionNodes * 8),
        DataType:  kernel_program.Float64,
        Alignment: kernel_program.CacheLineAlign,
    },
    {
        Name:      "rhs1",
        Size:      int64(totalSolutionNodes * 8),
        DataType:  kernel_program.Float64,
        Alignment: kernel_program.CacheLineAlign,
    },
    // ... rhs2, rhs3, rhs4
    
    // Derivative arrays for Dr, Ds, Dt operations
    {
        Name:      "ur",
        Size:      int64(totalSolutionNodes * 8),
        DataType:  kernel_program.Float64,
        Alignment: kernel_program.CacheLineAlign,
    },
    {
        Name:      "us", 
        Size:      int64(totalSolutionNodes * 8),
        DataType:  kernel_program.Float64,
        Alignment: kernel_program.CacheLineAlign,
    },
    {
        Name:      "ut",
        Size:      int64(totalSolutionNodes * 8),
        DataType:  kernel_program.Float64,
        Alignment: kernel_program.CacheLineAlign,
    },
    
    // Flux correction array
    {
        Name:      "flux_correction",
        Size:      int64(totalFaceNodes * 8),
        DataType:  kernel_program.Float64,
        Alignment: kernel_program.CacheLineAlign,
    },
    
    // Face geometry arrays
    {
        Name:      "nx",
        Size:      int64(totalFaceNodes * 8),
        DataType:  kernel_program.Float64,
        Alignment: kernel_program.CacheLineAlign,
    },
    {
        Name:      "ny",
        Size:      int64(totalFaceNodes * 8),
        DataType:  kernel_program.Float64,
        Alignment: kernel_program.CacheLineAlign,
    },
    {
        Name:      "nz",
        Size:      int64(totalFaceNodes * 8),
        DataType:  kernel_program.Float64,
        Alignment: kernel_program.CacheLineAlign,
    },
    {
        Name:      "Fscale",
        Size:      int64(totalFaceNodes * 8),
        DataType:  kernel_program.Float64,
        Alignment: kernel_program.CacheLineAlign,
    },
    // P_indices for neighbor connectivity
    {
        Name:      "P_indices",
        Size:      int64(totalFaceNodes * 8), // int64 size
        DataType:  kernel_program.INT64,
        Alignment: kernel_program.CacheLineAlign,
    },
}

// Allocate all arrays
err = kp.AllocateArrays(arrays)
if err != nil {
    log.Fatalf("Failed to allocate arrays: %v", err)
}
```

### Step 6: Initialize Arrays from Split Mesh

```go
// Helper to pack partition data
func packPartitionData(elements []*tetelement.Element3D, field func(el *tetelement.Element3D) []float64) []float64 {
    var packed []float64
    for _, partEl := range elements {
        packed = append(packed, field(partEl)...)
    }
    return packed
}

// Initialize solution with test function
initialU := make([]float64, totalSolutionNodes)
idx := 0
for _, partEl := range workingElements {
    for k := 0; k < partEl.K; k++ {
        for n := 0; n < partEl.Np; n++ {
            x := partEl.X.DataP[k*partEl.Np + n]
            y := partEl.Y.DataP[k*partEl.Np + n]
            z := partEl.Z.DataP[k*partEl.Np + n]
            // Gaussian pulse initial condition
            initialU[idx] = math.Exp(-10 * (x*x + y*y + z*z))
            idx++
        }
    }
}

// Copy arrays to device
kp.CopyArrayToDevice("u", initialU)
kp.CopyArrayToDevice("Rx", packPartitionData(workingElements, func(e *tetelement.Element3D) []float64 { return e.Rx.DataP }))
kp.CopyArrayToDevice("Ry", packPartitionData(workingElements, func(e *tetelement.Element3D) []float64 { return e.Ry.DataP }))
// ... repeat for all geometric factors and face arrays
```

### Step 7: Create Face Exchange Buffers

```go
// Structure to manage face exchanges between partitions
type FaceBuffer struct {
    sendBuffer []float64
    recvBuffer []float64
    sendMap    []int  // indices to extract from M
    recvMap    []int  // indices to place in P_buffer
    neighbor   int    // partition ID
}

// Build face exchange maps from EToE across partitions
faceBuffers := buildFaceExchangeBuffers(workingElements)
```

### Step 8: Build Computational Kernels

The matrix multiplication macros take (IN, OUT, K_VAL) parameters. Matrix dimensions are automatically inferred from the static matrix dimensions. The @inner loop is contained within the macro. MATMUL_ADD_* variants accumulate to output for flux corrections.

```go
// Volume RHS kernel for scalar Burger's equation
volumeKernelCode := `
#define NP ` + fmt.Sprintf("%d", Np) + `

@kernel void volumeRHS3D(
    const int_t* K,
    const real_t* u_global, const int_t* u_offsets,
    real_t* ur_global, const int_t* ur_offsets,
    real_t* us_global, const int_t* us_offsets,
    real_t* ut_global, const int_t* ut_offsets,
    const real_t* Rx_global, const int_t* Rx_offsets,
    const real_t* Ry_global, const int_t* Ry_offsets,
    const real_t* Rz_global, const int_t* Rz_offsets,
    const real_t* Sx_global, const int_t* Sx_offsets,
    const real_t* Sy_global, const int_t* Sy_offsets,
    const real_t* Sz_global, const int_t* Sz_offsets,
    const real_t* Tx_global, const int_t* Tx_offsets,
    const real_t* Ty_global, const int_t* Ty_offsets,
    const real_t* Tz_global, const int_t* Tz_offsets,
    real_t* rhs_global, const int_t* rhs_offsets
) {
    for (int part = 0; part < NPART; ++part; @outer) {
        const real_t* u = u_PART(part);
        real_t* ur = ur_PART(part);
        real_t* us = us_PART(part);
        real_t* ut = ut_PART(part);
        const real_t* Rx = Rx_PART(part);
        const real_t* Ry = Ry_PART(part);
        const real_t* Rz = Rz_PART(part);
        const real_t* Sx = Sx_PART(part);
        const real_t* Sy = Sy_PART(part);
        const real_t* Sz = Sz_PART(part);
        const real_t* Tx = Tx_PART(part);
        const real_t* Ty = Ty_PART(part);
        const real_t* Tz = Tz_PART(part);
        real_t* rhs = rhs_PART(part);
        
        int k_part = K[part];
        
        // Apply differentiation matrices
        MATMUL_Dr(u, ur, k_part);
        MATMUL_Ds(u, us, k_part);
        MATMUL_Dt(u, ut, k_part);
        
        // Compute flux divergence
        for (int elem = 0; elem < KpartMax; ++elem; @inner) {
            if (elem < k_part) {
                for (int i = 0; i < NP; ++i) {
                    int gid = elem * NP + i;
                    
                    // Get solution value
                    real_t u_local = u[gid];
                    
                    // Compute physical derivatives
                    real_t ux = Rx[gid]*ur[gid] + Sx[gid]*us[gid] + Tx[gid]*ut[gid];
                    real_t uy = Ry[gid]*ur[gid] + Sy[gid]*us[gid] + Ty[gid]*ut[gid];
                    real_t uz = Rz[gid]*ur[gid] + Sz[gid]*us[gid] + Tz[gid]*ut[gid];
                    
                    // Flux derivatives: ∇·F where F = (u²/2, u²/2, u²/2)
                    // dF/dx = u*ux, dG/dy = u*uy, dH/dz = u*uz
                    real_t dFdx = u_local * ux;
                    real_t dGdy = u_local * uy;
                    real_t dHdz = u_local * uz;
                    
                    // Volume contribution: -∇·F
                    rhs[gid] = -(dFdx + dGdy + dHdz);
                }
            }
        }
    }
}
`

// Extract face values kernel
extractFaceKernelCode := `
#define NP ` + fmt.Sprintf("%d", Np) + `
#define NFP ` + fmt.Sprintf("%d", Nfp) + `

@kernel void extractFaceValues(
    const int_t* K,
    const real_t* u_global, const int_t* u_offsets,
    real_t* M_global, const int_t* M_offsets,
    const int_t* VmapM
) {
    for (int part = 0; part < NPART; ++part; @outer) {
        const real_t* u = u_PART(part);
        real_t* M = M_PART(part);
        int k_part = K[part];
        
        // Extract values at face nodes
        for (int elem = 0; elem < KpartMax; ++elem; @inner) {
            if (elem < k_part) {
                for (int f = 0; f < 4; ++f) {
                    for (int i = 0; i < NFP; ++i) {
                        int fid = elem * NFP * 4 + f * NFP + i;
                        int vid = VmapM[f * NFP + i] + elem * NP;
                        M[fid] = u[vid];
                    }
                }
            }
        }
    }
}
`

// Flux correction kernel
fluxCorrectionKernelCode := `
#define NP ` + fmt.Sprintf("%d", Np) + `
#define NFP ` + fmt.Sprintf("%d", Nfp) + `

@kernel void fluxCorrection3D(
    const int_t* K,
    const real_t* M_global, const int_t* M_offsets,
    const real_t* nx_global, const int_t* nx_offsets,
    const real_t* ny_global, const int_t* ny_offsets,
    const real_t* nz_global, const int_t* nz_offsets,
    const real_t* Fscale_global, const int_t* Fscale_offsets,
    const int_t* P_indices_global, const int_t* P_indices_offsets,
    real_t* flux_correction_global, const int_t* flux_correction_offsets,
    real_t* rhs_global, const int_t* rhs_offsets
) {
    for (int part = 0; part < NPART; ++part; @outer) {
        const real_t* M = M_PART(part);
        const real_t* nx = nx_PART(part);
        const real_t* ny = ny_PART(part);
        const real_t* nz = nz_PART(part);
        const real_t* Fscale = Fscale_PART(part);
        const int_t* P_indices = P_indices_PART(part);
        real_t* flux_correction = flux_correction_PART(part);
        real_t* rhs = rhs_PART(part);
        
        int k_part = K[part];
        
        // Compute flux corrections at faces
        for (int elem = 0; elem < KpartMax; ++elem; @inner) {
            if (elem < k_part) {
                for (int f = 0; f < 4; ++f) {
                    for (int i = 0; i < NFP; ++i) {
                        int fid = elem * NFP * 4 + f * NFP + i;
                        real_t uM = M[fid];
                        
                        // Get uP from P buffer location
                        int p_idx = P_indices[fid];
                        real_t uP = (p_idx >= 0) ? M[p_idx] : uM; // Boundary: uP = uM
                        
                        // Local flux (conservation form)
                        real_t FM = 0.5 * uM * uM;
                        real_t GM = FM;  // Same for all directions in Burger's
                        real_t HM = FM;
                        
                        // Numerical flux (Lax-Friedrichs)
                        real_t alpha = fmax(fabs(uM), fabs(uP));
                        real_t Fnum = 0.5 * (0.5*(uM*uM + uP*uP) - alpha*(uP - uM));
                        real_t Gnum = Fnum;
                        real_t Hnum = Fnum;
                        
                        // Normal flux: (F* - F)·n̂
                        real_t flux_jump = (Fnum - FM)*nx[fid] + 
                                          (Gnum - GM)*ny[fid] + 
                                          (Hnum - HM)*nz[fid];
                        
                        // Scale by face Jacobian
                        flux_correction[elem * NFP * 4 + f * NFP + i] = flux_jump * Fscale[fid];
                    }
                }
            }
        }
        
        // Apply LIFT matrix to get volume contribution
        // LIFT is [Np x Nfp*4], flux_correction is [Nfp*4 x K]
        // Result accumulates to rhs [Np x K]
        MATMUL_ADD_LIFT(flux_correction, rhs, k_part);
    }
}
`

// Build kernels
kp.BuildKernel("volumeRHS3D", volumeKernelCode)
kp.BuildKernel("extractFaceValues", extractFaceKernelCode)
kp.BuildKernel("fluxCorrection3D", fluxCorrectionKernelCode)
```

### Step 9: Time Integration Loop

```go
// RK5SSP4 coefficients
c := []float64{0, 0.391752226571890, 0.555629506348765, 0.682007177596368, 1.0}
a := [][]float64{
    {0},
    {0.391752226571890},
    {0.216881479655892, 0.338748077291873},
    {0.082637967089070, 0.129460862868150, 0.469908347471148},
    {0.066740164024413, 0.104615396254984, 0.379215311156271, 0.449429027564318},
}

// Copy VmapM to device once (static data)
// VmapM maps face nodes to volume nodes
vmapMData := make([]int32, totalFaceNodes)
idx = 0
for _, partEl := range workingElements {
    for k := 0; k < partEl.K; k++ {
        for f := 0; f < 4; f++ {
            for i := 0; i < partEl.Nfp; i++ {
                // VmapM maps face node to volume node
                vmapMData[idx] = int32(partEl.VmapM[k*4*partEl.Nfp + f*partEl.Nfp + i])
                idx++
            }
        }
    }
}
vmapMDevice := device.Malloc(int64(len(vmapMData)*4), unsafe.Pointer(&vmapMData[0]))
defer vmapMDevice.Free()

// Main time loop
for step := 0; step < numSteps; step++ {
    // RK5SSP4 stages
    for stage := 0; stage < 5; stage++ {
        // 1. Extract face values to M buffer
        err = kp.RunKernelWithExtra("extractFaceValues", 
            []interface{}{vmapMDevice}, 
            "u", "M")
        
        // 2. Exchange partition boundaries
        exchangePartitionBoundaries(kp, faceBuffers)
        
        // 3. Compute volume RHS (flux divergence)
        rhsName := fmt.Sprintf("rhs%d", stage)
        err = kp.RunKernel("volumeRHS3D", 
            "u", "ur", "us", "ut", "Rx", "Ry", "Rz", "Sx", "Sy", "Sz", "Tx", "Ty", "Tz", rhsName)
        
        // 4. Compute flux corrections and apply LIFT
        err = kp.RunKernel("fluxCorrection3D",
            "M", "nx", "ny", "nz", "Fscale", "P_indices", "BC_types",
            "flux_correction", rhsName)
        
        // 5. RK update (stage-specific linear combination)
        if stage == 0 {
            // u = u + dt * a[stage][0] * rhs0
            kp.RunKernel("rkUpdate", "u", "rhs0", dt*a[stage][0])
        } else {
            // u = weighted sum of previous stages
            // Implementation depends on stage
        }
    }
    
    // Monitor solution
    if step % 100 == 0 {
        u_host, _ := kp.CopyArrayToHost("u")
        fmt.Printf("Step %d: max|u| = %e\n", step, maxAbs(u_host))
    }
}
```

### Step 10: Post-Processing

```go
// Copy solution back to host
finalU, err := kp.CopyArrayToHost("u")
if err != nil {
    log.Fatalf("Failed to copy solution: %v", err)
}

// Unpack from partition format to element format
// Write VTK output using existing Element3D methods
```

## Key Implementation Notes

### 1. Memory Layout

- Arrays are allocated as contiguous blocks with automatic offset management
- KernelProgram handles partition access through generated macros
- Derivative arrays (ur, us, ut) allocated as device memory to avoid stack issues
- Alignment ensures optimal memory access patterns

### 2. Kernel Design

- All kernels follow OCCA’s @outer/@inner pattern
- Matrix operations use generated macros containing @inner loops
- MATMUL_* macros take only (IN, OUT, K_VAL) - matrix dimensions are embedded
- MATMUL_ADD_* variants accumulate to output (needed for flux corrections)
- Flux correction array allocated as device memory to avoid stack issues

### 3. Partition Parallelism

- Each partition processes independently in @outer loop
- Face exchanges occur between RHS evaluations
- Variable partition sizes handled through K array

### 4. Numerical Methods

- Conservation form with nonlinear flux (F = u²/2)
- Lax-Friedrichs numerical flux for stability
- Wall BC: zero normal flux by removing **F**M·**n̂** component
- LIFT matrix maps face corrections to volume contributions

### 5. Performance Considerations

- Static matrices embedded in kernel code
- Memory alignment for cache efficiency
- Coalesced access patterns in kernels
- CUDA thread limit checked during initialization

## Testing Strategy

Following the Unit Testing Principles:

1. **Start with fundamentals**: Test single element, then single partition
1. **Build systematically**: Add partitions incrementally
1. **Specific properties**: Verify conservation, convergence rates
1. **Detailed coverage**: Test boundary conditions, parallel exchange

### Test Cases

1. **Single Element Tests**
- Constant solution preservation
- Linear advection accuracy (before shock formation)
- Conservation properties
1. **Multi-Partition Tests**
- Partition boundary continuity
- Load balancing verification
- Parallel efficiency
1. **Convergence Tests**
- Smooth solution: O(N+1) convergence
- Discontinuous solution: shock capturing
1. **Benchmark Problems**
- 3D advection of Gaussian pulse
- Shock formation from smooth initial data
- Interaction of multiple shocks
- Rarefaction wave propagation

## Error Handling

- Validate mesh connectivity before processing
- Check CUDA limits for partition sizes
- Verify geometric factor positivity (J > 0)
- Monitor CFL condition during time stepping
- Handle array allocation failures gracefully

## Extensions to Navier-Stokes

This Burger’s equation solver demonstrates key capabilities needed for Navier-Stokes:

1. **Flux Formulation**: The conservation form with flux derivatives directly extends to the convective terms in Navier-Stokes
1. **Face Corrections**: The flux correction approach handles discontinuous solutions needed for compressible flows
1. **Partition Parallel**: The infrastructure supports the computational intensity of 3D Navier-Stokes
1. **Geometric Flexibility**: Unstructured tetrahedral meshes handle complex geometries

The main additions for Navier-Stokes will be:

- Vector solutions (ρ, ρu, ρv, ρw, E)
- Viscous flux terms requiring gradient computation
- Boundary condition variety (wall, inflow, outflow)
- Shock capturing and limiting strategies

## Boundary Condition Implementation

The solver supports wall and farfield boundary conditions through the flux computation. Boundary conditions are implemented in the flux correction kernel by modifying how the neighbor value uP is determined.

### Boundary Condition Mapping

```go
// Build BC mapping during initialization
// BCMaps contains face-to-BC-type mapping from mesh
bcTypes := make([]int32, totalFaceNodes)
idx := 0
for _, partEl := range workingElements {
    for k := 0; k < partEl.K; k++ {
        for f := 0; f < 4; f++ {
            for i := 0; i < partEl.Nfp; i++ {
                faceIdx := k*4*partEl.Nfp + f*partEl.Nfp + i
                
                // Check if this is a boundary face
                if partEl.EToE[k][f] == k { // Self-connected = boundary
                    // Get BC type from mesh
                    bcTypes[idx] = int32(partEl.GetFaceBCType(k, f))
                } else {
                    bcTypes[idx] = BC_NONE // Internal face
                }
                idx++
            }
        }
    }
}

// Copy to device
bcTypesDevice := device.Malloc(int64(len(bcTypes)*4), unsafe.Pointer(&bcTypes[0]))
defer bcTypesDevice.Free()
```

### Updated Flux Correction Kernel with BCs

```c
// Add BC type constants
#define BC_NONE 0
#define BC_WALL 1
#define BC_FARFIELD 2

@kernel void fluxCorrection3D(
    const int_t* K,
    const real_t* M_global, const int_t* M_offsets,
    const real_t* nx_global, const int_t* nx_offsets,
    const real_t* ny_global, const int_t* ny_offsets,
    const real_t* nz_global, const int_t* nz_offsets,
    const real_t* Fscale_global, const int_t* Fscale_offsets,
    const int_t* P_indices_global, const int_t* P_indices_offsets,
    const int_t* BC_types_global, const int_t* BC_types_offsets,
    real_t* flux_correction_global, const int_t* flux_correction_offsets,
    real_t* rhs_global, const int_t* rhs_offsets
) {
    for (int part = 0; part < NPART; ++part; @outer) {
        const real_t* M = M_PART(part);
        const real_t* nx = nx_PART(part);
        const real_t* ny = ny_PART(part);
        const real_t* nz = nz_PART(part);
        const real_t* Fscale = Fscale_PART(part);
        const int_t* P_indices = P_indices_PART(part);
        const int_t* BC_types = BC_types_PART(part);
        real_t* flux_correction = flux_correction_PART(part);
        real_t* rhs = rhs_PART(part);
        
        int k_part = K[part];
        
        // Compute flux corrections at faces
        for (int elem = 0; elem < KpartMax; ++elem; @inner) {
            if (elem < k_part) {
                for (int f = 0; f < 4; ++f) {
                    for (int i = 0; i < NFP; ++i) {
                        int fid = elem * NFP * 4 + f * NFP + i;
                        real_t uM = M[fid];
                        real_t uP;
                        
                        // Get boundary condition type
                        int bc_type = BC_types[fid];
                        
                        if (bc_type == BC_NONE) {
                            // Internal face - get neighbor value
                            int p_idx = P_indices[fid];
                            uP = (p_idx >= 0) ? M[p_idx] : uM;
                        } else if (bc_type == BC_WALL) {
                            // Wall BC: zero normal flux condition
                            // For zero flux, the numerical flux f* must equal fM
                            // so that (f* - fM)·n̂ = 0
                            // This is achieved by setting uP = uM (not -uM!)
                            uP = uM;
                        } else if (bc_type == BC_FARFIELD) {
                            // Farfield BC: Riemann invariant
                            // For outflow (u·n > 0): extrapolate from interior
                            // For inflow (u·n < 0): use farfield value
                            real_t un = uM * 1.0; // Simplified for scalar
                            if (un > 0) {
                                uP = uM; // Outflow: extrapolate
                            } else {
                                uP = 0.0; // Inflow: farfield value (customize as needed)
                            }
                        } else {
                            // Default: treat as wall
                            uP = uM;
                        }
                        
                        // Local flux (conservation form)
                        real_t FM = 0.5 * uM * uM;
                        real_t FP = 0.5 * uP * uP;
                        
                        // Numerical flux (Lax-Friedrichs)
                        real_t alpha = fmax(fabs(uM), fabs(uP));
                        real_t Fnum = 0.5 * (FM + FP - alpha*(uP - uM));
                        
                        // For scalar Burger's, flux is same in all directions
                        real_t Gnum = Fnum;
                        real_t Hnum = Fnum;
                        
                        // Normal flux: (F* - F)·n̂
                        real_t flux_jump = (Fnum - FM)*nx[fid] + 
                                          (Gnum - FM)*ny[fid] + 
                                          (Hnum - FM)*nz[fid];
                        
                        // Scale by face Jacobian
                        flux_correction[elem * NFP * 4 + f * NFP + i] = flux_jump * Fscale[fid];
                    }
                }
            }
        }
        
        // Apply LIFT matrix to get volume contribution
        MATMUL_ADD_LIFT(flux_correction, rhs, k_part);
    }
}
```

### Boundary Condition Notes

1. **Wall BC**: For the scalar Burger’s equation, the wall boundary condition enforces zero normal flux: **F**·**n̂** = 0. Since **F** = (u²/2, u²/2, u²/2), the normal flux is (u²/2)(nx + ny + nz). The flux correction at the wall is simply -**F**M·**n̂**, which removes the normal flux component.
1. **Farfield BC**: Uses characteristic analysis based on the wave speed:
- For outflow (u > 0): Information travels out, so we extrapolate from interior
- For inflow (u < 0): Information comes from outside, so we impose farfield values
1. **Extension to Navier-Stokes**:
- This scalar formulation tests the DG infrastructure and shock-capturing capabilities
- The flux correction approach extends directly to systems of equations
- Euler equations add pressure terms to momentum fluxes
- Navier-Stokes adds viscous flux terms requiring gradient computation

## References

- Hesthaven & Warburton, “Nodal Discontinuous Galerkin Methods”
- Spiteri & Ruuth, “Optimal Strong-Stability-Preserving Runge-Kutta Methods”
- GOCCA documentation for kernel programming patterns
- KernelProgram User’s Guide for array management

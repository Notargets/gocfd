# 3D Inviscid Burger’s Equation Solver Implementation Guide

## Overview

This document describes the implementation of a 3D inviscid Burger’s equation solver using discontinuous Galerkin (DG) methods in conservation form. The solver combines `kernel_program.go` from gocca with the tetrahedra package from gocfd to create a partition-parallel solver supporting both OpenMP and CUDA execution.

The conservation form approach with flux corrections serves as a testbed for the numerical methods and infrastructure that will be used in a subsequent Navier-Stokes equations solver.

## Mathematical Formulation

The 3D inviscid Burger’s equation in conservation form:

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
// Get partition data from split mesh
numPartitions := len(el.Split)
K := make([]int, numPartitions)
for i, partEl := range el.Split {
    K[i] = partEl.K
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
    
    // Normal vectors for surface integrals
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
    // Fscale for surface integrals
    {
        Name:      "Fscale",
        Size:      int64(totalFaceNodes * 8),
        DataType:  kernel_program.Float64,
        Alignment: kernel_program.CacheLineAlign,
    },
}

// Allocate all arrays
err = kp.AllocateArrays(arrays)
if err != nil {
    log.Fatal("Failed to allocate arrays:", err)
}
```

### Step 6: Copy Geometric Factors to Device

```go
// Helper function to copy partition data to device
func copyPartitionDataToDevice(kp *kernel_program.KernelProgram, 
                               arrayName string, 
                               partitions []*tetelement.Element3D,
                               getData func(*tetelement.Element3D) utils.Matrix) {
    mem := kp.GetMemory(arrayName)
    
    offset := 0
    for _, partEl := range partitions {
        matrix := getData(partEl)
        size := matrix.Len() * 8  // 8 bytes per float64
        
        // Copy matrix data
        mem.CopyFromWithOffset(unsafe.Pointer(matrix.DataP), int64(size), int64(offset))
        offset += size
    }
}

// Copy all geometric factors
copyPartitionDataToDevice(kp, "Rx", el.Split, func(e *tetelement.Element3D) utils.Matrix { return e.Rx })
copyPartitionDataToDevice(kp, "Ry", el.Split, func(e *tetelement.Element3D) utils.Matrix { return e.Ry })
// ... repeat for all geometric factors

// Copy normal vectors and Fscale
copyPartitionDataToDevice(kp, "nx", el.Split, func(e *tetelement.Element3D) utils.Matrix { return e.Nx })
copyPartitionDataToDevice(kp, "ny", el.Split, func(e *tetelement.Element3D) utils.Matrix { return e.Ny })
copyPartitionDataToDevice(kp, "nz", el.Split, func(e *tetelement.Element3D) utils.Matrix { return e.Nz })
copyPartitionDataToDevice(kp, "Fscale", el.Split, func(e *tetelement.Element3D) utils.Matrix { return e.Fscale })
```

### Step 7: Build Face Exchange P Indices

```go
// Build face buffers for each partition
faceBuffers := make([]*facebuffer.FaceBuffer, numPartitions)
totalPIndices := 0

for p, partEl := range el.Split {
    fb, err := facebuffer.BuildFaceBuffer(partEl)
    if err != nil {
        log.Fatal(err)
    }
    faceBuffers[p] = fb
    totalPIndices += len(fb.LocalPIndices)
}

// Allocate P indices array if needed
if totalPIndices > 0 {
    err = kp.AllocateArrays([]kernel_program.ArraySpec{
        {
            Name:      "P_indices",
            Size:      int64(totalPIndices * 4), // 4 bytes per int32
            DataType:  kernel_program.INT32,
            Alignment: kernel_program.NoAlignment,
        },
    })
    
    // Copy P indices to device
    mem := kp.GetMemory("P_indices")
    offset := 0
    for _, fb := range faceBuffers {
        if len(fb.LocalPIndices) > 0 {
            size := len(fb.LocalPIndices) * 4
            mem.CopyFromWithOffset(unsafe.Pointer(&fb.LocalPIndices[0]), int64(size), int64(offset))
            offset += size
        }
    }
}
```

### Step 8: Define Kernels

```go
// Volume RHS kernel - computes flux divergence
volumeKernelCode := `
#define NP ` + fmt.Sprintf("%d", Np) + `
#define NFP ` + fmt.Sprintf("%d", Nfp) + `

@kernel void volumeRHS3D(
    const int_t* K,
    const real_t* u_global, const int_t* u_offsets,
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
        // Get partition data pointers
        const real_t* u = u_PART(part);
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
        real_t ur[NP * KpartMax];
        real_t us[NP * KpartMax];
        real_t ut[NP * KpartMax];
        
        MATMUL_Dr(u, ur, k_part, KpartMax, NP);
        MATMUL_Ds(u, us, k_part, KpartMax, NP);
        MATMUL_Dt(u, ut, k_part, KpartMax, NP);
        
        // Compute flux derivatives
        for (int elem = 0; elem < KpartMax; ++elem; @inner) {
            if (elem < k_part) {
                for (int i = 0; i < NP; ++i) {
                    int gid = elem * NP + i;
                    
                    // Get solution value
                    real_t u_local = u[gid];
                    real_t flux = 0.5 * u_local * u_local;
                    
                    // Compute flux derivatives using chain rule
                    // dF/dx = dF/du * du/dx = u * ux
                    real_t ux = Rx[gid]*ur[gid] + Sx[gid]*us[gid] + Tx[gid]*ut[gid];
                    real_t uy = Ry[gid]*ur[gid] + Sy[gid]*us[gid] + Ty[gid]*ut[gid];
                    real_t uz = Rz[gid]*ur[gid] + Sz[gid]*us[gid] + Tz[gid]*ut[gid];
                    
                    real_t dFdx = u_local * ux;
                    real_t dGdy = u_local * uy;
                    real_t dHdz = u_local * uz;
                    
                    // Volume contribution: -div(F)
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
    real_t* rhs_global, const int_t* rhs_offsets
) {
    for (int part = 0; part < NPART; ++part; @outer) {
        const real_t* M = M_PART(part);
        const real_t* nx = nx_PART(part);
        const real_t* ny = ny_PART(part);
        const real_t* nz = nz_PART(part);
        const real_t* Fscale = Fscale_PART(part);
        const int_t* P_indices = P_indices_PART(part);
        real_t* rhs = rhs_PART(part);
        
        int k_part = K[part];
        
        // Compute flux corrections at faces
        for (int elem = 0; elem < KpartMax; ++elem; @inner) {
            if (elem < k_part) {
                real_t flux_correction[NFP * 4];
                
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
                        flux_correction[f * NFP + i] = flux_jump * Fscale[fid];
                    }
                }
                
                // Apply LIFT matrix to get volume contribution
                for (int i = 0; i < NP; ++i) {
                    real_t lift_contribution = 0.0;
                    for (int j = 0; j < NFP * 4; ++j) {
                        lift_contribution += LIFT[i][j] * flux_correction[j];
                    }
                    @atomic rhs[elem * NP + i] += lift_contribution;
                }
            }
        }
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
vmapMData := make([]int32, Nfp*4)
for f := 0; f < 4; f++ {
    for i := 0; i < Nfp; i++ {
        vmapMData[f*Nfp + i] = int32(el.Fmask[f][i])
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
            "u", "Rx", "Ry", "Rz", "Sx", "Sy", "Sz", "Tx", "Ty", "Tz", rhsName)
        
        // 4. Compute flux corrections and apply LIFT
        err = kp.RunKernel("fluxCorrection3D",
            "M", "nx", "ny", "nz", "Fscale", "P_indices", rhsName)
        
        // 5. Update solution (implement as kernel)
        updateSolution(kp, stage, dt, a, c)
    }
    
    // Monitor solution
    if step % 100 == 0 {
        maxU := computeMaximum(kp, "u")
        fmt.Printf("Step %d: max(u) = %g\n", step, maxU)
    }
}
```

### Step 10: Helper Functions

```go
// Exchange partition boundary data
func exchangePartitionBoundaries(kp *kernel_program.KernelProgram, 
                                faceBuffers []*facebuffer.FaceBuffer) {
    // For each partition, copy remote face data
    for p, fb := range faceBuffers {
        for remotePartID, indices := range fb.RemoteSendIndices {
            // Copy from partition p to remotePartID
            // This requires MPI or similar for distributed memory
            // For shared memory, direct copy between M buffers
        }
    }
}

// Compute maximum value (for monitoring)
func computeMaximum(kp *kernel_program.KernelProgram, arrayName string) float64 {
    data, _ := kernel_program.CopyArrayToHost[float64](kp, arrayName)
    maxVal := data[0]
    for _, v := range data {
        if v > maxVal {
            maxVal = v
        }
    }
    return maxVal
}

// Update solution kernel
func buildUpdateKernel(kp *kernel_program.KernelProgram, stage int, a [][]float64) {
    code := fmt.Sprintf(`
    @kernel void updateStage%d(
        const int_t* K,
        real_t* u_global, const int_t* u_offsets,
        const real_t* u0_global, const int_t* u0_offsets,
    `, stage)
    
    // Add RHS arrays used in this stage
    for s := 0; s <= stage; s++ {
        code += fmt.Sprintf(`    const real_t* rhs%d_global, const int_t* rhs%d_offsets,
    `, s, s)
    }
    
    code += fmt.Sprintf(`    const real_t dt
    ) {
        for (int part = 0; part < NPART; ++part; @outer) {
            real_t* u = u_PART(part);
            const real_t* u0 = u0_PART(part);
            `)
    
    // Add RHS pointers
    for s := 0; s <= stage; s++ {
        code += fmt.Sprintf(`const real_t* rhs%d = rhs%d_PART(part);
            `, s, s)
    }
    
    code += `
            int k_part = K[part];
            
            for (int i = 0; i < k_part * NP; ++i; @inner) {
                u[i] = u0[i]`
    
    // Add RHS contributions
    for s := 0; s <= stage; s++ {
        if a[stage][s] > 0 {
            code += fmt.Sprintf(` + dt * %g * rhs%d[i]`, a[stage][s], s)
        }
    }
    
    code += `;
            }
        }
    }`
    
    kp.BuildKernel(fmt.Sprintf("updateStage%d", stage), code)
}
```

## Key Implementation Details

### 1. Memory Layout

- Arrays are allocated as contiguous blocks with automatic offset management
- KernelProgram handles partition access through generated macros
- Alignment ensures optimal memory access patterns

### 2. Kernel Design

- All kernels follow OCCA’s @outer/@inner pattern
- Matrix operations use generated macros containing @inner loops
- Atomic operations handle race conditions in LIFT application

### 3. Partition Parallelism

- Each partition processes independently in @outer loop
- Face exchanges occur between RHS evaluations
- Variable partition sizes handled through K array

### 4. Numerical Methods

- Conservation form with flux formulation (F = u²/2)
- Lax-Friedrichs numerical flux for stability
- Flux corrections computed as (F* - F)·n̂
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
- Linear advection accuracy
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

## References

- Hesthaven & Warburton, “Nodal Discontinuous Galerkin Methods”
- Spiteri & Ruuth, “Optimal Strong-Stability-Preserving Runge-Kutta Methods”
- GOCCA documentation for kernel programming patterns
- KernelProgram User’s Guide for array management

# 3D Inviscid Burger's Equation Solver Implementation Guide

## Overview

This document describes the implementation of a 3D vector inviscid Burger's equation solver using discontinuous Galerkin (DG) methods in conservation form. The solver combines `kernel_program.go` from gocca with the tetrahedra package from gocfd to create a partition-parallel solver supporting both OpenMP and CUDA execution.

The conservation form approach with flux corrections serves as a testbed for the numerical methods and infrastructure that will be used in a subsequent Navier-Stokes equations solver.

## Mathematical Formulation

The 3D vector inviscid Burger's equation in non-conservative form:

```
∂u/∂t + (u·∇)u = 0

where u = (u, v, w) is the velocity vector
```

In conservation form, this becomes three coupled equations:

```
∂u/∂t + ∂(u²)/∂x + ∂(uv)/∂y + ∂(uw)/∂z = 0
∂v/∂t + ∂(uv)/∂x + ∂(v²)/∂y + ∂(vw)/∂z = 0
∂w/∂t + ∂(uw)/∂x + ∂(vw)/∂y + ∂(w²)/∂z = 0
```

The flux tensor is:
```
F = [u², uv, uw]ᵀ    (x-direction fluxes)
G = [uv, v², vw]ᵀ    (y-direction fluxes)
H = [uw, vw, w²]ᵀ    (z-direction fluxes)
```

In the DG formulation:

```
∂q/∂t = -∇·F^v + LIFT * (F* - F)·n̂

where:
- q = (u, v, w) is the solution vector
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
    totalSolutionNodes += Np * k * 3  // 3 components (u, v, w)
    totalFaceNodes += Nfp * 4 * k * 3  // 4 faces per tet, 3 components
    totalGeometricFactors += Np * k
}

// Define array specifications
arrays := []kernel_program.ArraySpec{
    // Solution arrays - 3 components
    {
        Name:      "u",
        Size:      int64(totalSolutionNodes * 8), // 8 bytes per float64
        DataType:  kernel_program.Float64,
        Alignment: kernel_program.CacheLineAlign,
    },
    {
        Name:      "v",
        Size:      int64(totalSolutionNodes/3 * 8), // Single component size
        DataType:  kernel_program.Float64,
        Alignment: kernel_program.CacheLineAlign,
    },
    {
        Name:      "w",
        Size:      int64(totalSolutionNodes/3 * 8),
        DataType:  kernel_program.Float64,
        Alignment: kernel_program.CacheLineAlign,
    },
    // M buffer for face data - 3 components
    {
        Name:      "Mu",
        Size:      int64(totalFaceNodes/3 * 8),
        DataType:  kernel_program.Float64,
        Alignment: kernel_program.CacheLineAlign,
    },
    {
        Name:      "Mv",
        Size:      int64(totalFaceNodes/3 * 8),
        DataType:  kernel_program.Float64,
        Alignment: kernel_program.CacheLineAlign,
    },
    {
        Name:      "Mw",
        Size:      int64(totalFaceNodes/3 * 8),
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
    
    // RHS arrays for RK5SSP4 (5 stages) - 3 components each
    {
        Name:      "rhsu0",
        Size:      int64(totalSolutionNodes/3 * 8),
        DataType:  kernel_program.Float64,
        Alignment: kernel_program.CacheLineAlign,
    },
    {
        Name:      "rhsv0",
        Size:      int64(totalSolutionNodes/3 * 8),
        DataType:  kernel_program.Float64,
        Alignment: kernel_program.CacheLineAlign,
    },
    {
        Name:      "rhsw0",
        Size:      int64(totalSolutionNodes/3 * 8),
        DataType:  kernel_program.Float64,
        Alignment: kernel_program.CacheLineAlign,
    },
    // ... repeat for stages 1-4
    
    // Derivative arrays for Dr, Ds, Dt operations
    {
        Name:      "ur",
        Size:      int64(totalSolutionNodes/3 * 8),
        DataType:  kernel_program.Float64,
        Alignment: kernel_program.CacheLineAlign,
    },
    {
        Name:      "us", 
        Size:      int64(totalSolutionNodes/3 * 8),
        DataType:  kernel_program.Float64,
        Alignment: kernel_program.CacheLineAlign,
    },
    {
        Name:      "ut",
        Size:      int64(totalSolutionNodes/3 * 8),
        DataType:  kernel_program.Float64,
        Alignment: kernel_program.CacheLineAlign,
    },
    // ... repeat for vr, vs, vt, wr, ws, wt
    
    // Flux correction arrays
    {
        Name:      "flux_u",
        Size:      int64(totalFaceNodes/3 * 8),
        DataType:  kernel_program.Float64,
        Alignment: kernel_program.CacheLineAlign,
    },
    {
        Name:      "flux_v",
        Size:      int64(totalFaceNodes/3 * 8),
        DataType:  kernel_program.Float64,
        Alignment: kernel_program.CacheLineAlign,
    },
    {
        Name:      "flux_w",
        Size:      int64(totalFaceNodes/3 * 8),
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
initialU := make([]float64, totalSolutionNodes/3)
initialV := make([]float64, totalSolutionNodes/3)
initialW := make([]float64, totalSolutionNodes/3)
idx := 0
for _, partEl := range workingElements {
    for k := 0; k < partEl.K; k++ {
        for n := 0; n < partEl.Np; n++ {
            x := partEl.X.DataP[k*partEl.Np + n]
            y := partEl.Y.DataP[k*partEl.Np + n]
            z := partEl.Z.DataP[k*partEl.Np + n]
            // Example: rotating vortex initial condition
            r := math.Sqrt(x*x + y*y)
            if r > 0 {
                initialU[idx] = -y/r * math.Exp(-10 * (r*r + z*z))
                initialV[idx] = x/r * math.Exp(-10 * (r*r + z*z))
            } else {
                initialU[idx] = 0
                initialV[idx] = 0
            }
            initialW[idx] = 0
            idx++
        }
    }
}

// Copy arrays to device
kp.CopyArrayToDevice("u", initialU)
kp.CopyArrayToDevice("v", initialV)
kp.CopyArrayToDevice("w", initialW)
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
// Volume RHS kernel for vector Burger's equation
volumeKernelCode := `
#define NP ` + fmt.Sprintf("%d", Np) + `

@kernel void volumeRHS3D(
    const int_t* K,
    const real_t* u_global, const int_t* u_offsets,
    const real_t* v_global, const int_t* v_offsets,
    const real_t* w_global, const int_t* w_offsets,
    real_t* ur_global, const int_t* ur_offsets,
    real_t* us_global, const int_t* us_offsets,
    real_t* ut_global, const int_t* ut_offsets,
    real_t* vr_global, const int_t* vr_offsets,
    real_t* vs_global, const int_t* vs_offsets,
    real_t* vt_global, const int_t* vt_offsets,
    real_t* wr_global, const int_t* wr_offsets,
    real_t* ws_global, const int_t* ws_offsets,
    real_t* wt_global, const int_t* wt_offsets,
    const real_t* Rx_global, const int_t* Rx_offsets,
    const real_t* Ry_global, const int_t* Ry_offsets,
    const real_t* Rz_global, const int_t* Rz_offsets,
    const real_t* Sx_global, const int_t* Sx_offsets,
    const real_t* Sy_global, const int_t* Sy_offsets,
    const real_t* Sz_global, const int_t* Sz_offsets,
    const real_t* Tx_global, const int_t* Tx_offsets,
    const real_t* Ty_global, const int_t* Ty_offsets,
    const real_t* Tz_global, const int_t* Tz_offsets,
    real_t* rhsu_global, const int_t* rhsu_offsets,
    real_t* rhsv_global, const int_t* rhsv_offsets,
    real_t* rhsw_global, const int_t* rhsw_offsets
) {
    for (int part = 0; part < NPART; ++part; @outer) {
        const real_t* u = u_PART(part);
        const real_t* v = v_PART(part);
        const real_t* w = w_PART(part);
        real_t* ur = ur_PART(part);
        real_t* us = us_PART(part);
        real_t* ut = ut_PART(part);
        real_t* vr = vr_PART(part);
        real_t* vs = vs_PART(part);
        real_t* vt = vt_PART(part);
        real_t* wr = wr_PART(part);
        real_t* ws = ws_PART(part);
        real_t* wt = wt_PART(part);
        const real_t* Rx = Rx_PART(part);
        const real_t* Ry = Ry_PART(part);
        const real_t* Rz = Rz_PART(part);
        const real_t* Sx = Sx_PART(part);
        const real_t* Sy = Sy_PART(part);
        const real_t* Sz = Sz_PART(part);
        const real_t* Tx = Tx_PART(part);
        const real_t* Ty = Ty_PART(part);
        const real_t* Tz = Tz_PART(part);
        real_t* rhsu = rhsu_PART(part);
        real_t* rhsv = rhsv_PART(part);
        real_t* rhsw = rhsw_PART(part);
        
        int k_part = K[part];
        
        // Apply differentiation matrices
        MATMUL_Dr(u, ur, k_part);
        MATMUL_Ds(u, us, k_part);
        MATMUL_Dt(u, ut, k_part);
        
        MATMUL_Dr(v, vr, k_part);
        MATMUL_Ds(v, vs, k_part);
        MATMUL_Dt(v, vt, k_part);
        
        MATMUL_Dr(w, wr, k_part);
        MATMUL_Ds(w, ws, k_part);
        MATMUL_Dt(w, wt, k_part);
        
        // Compute flux divergence for vector Burger's equation
        for (int elem = 0; elem < KpartMax; ++elem; @inner) {
            if (elem < k_part) {
                for (int i = 0; i < NP; ++i) {
                    int gid = elem * NP + i;
                    
                    // Get velocity components
                    real_t u_local = u[gid];
                    real_t v_local = v[gid];
                    real_t w_local = w[gid];
                    
                    // Compute physical derivatives
                    real_t ux = Rx[gid]*ur[gid] + Sx[gid]*us[gid] + Tx[gid]*ut[gid];
                    real_t uy = Ry[gid]*ur[gid] + Sy[gid]*us[gid] + Ty[gid]*ut[gid];
                    real_t uz = Rz[gid]*ur[gid] + Sz[gid]*us[gid] + Tz[gid]*ut[gid];
                    
                    real_t vx = Rx[gid]*vr[gid] + Sx[gid]*vs[gid] + Tx[gid]*vt[gid];
                    real_t vy = Ry[gid]*vr[gid] + Sy[gid]*vs[gid] + Ty[gid]*vt[gid];
                    real_t vz = Rz[gid]*vr[gid] + Sz[gid]*vs[gid] + Tz[gid]*vt[gid];
                    
                    real_t wx = Rx[gid]*wr[gid] + Sx[gid]*ws[gid] + Tx[gid]*wt[gid];
                    real_t wy = Ry[gid]*wr[gid] + Sy[gid]*ws[gid] + Ty[gid]*wt[gid];
                    real_t wz = Rz[gid]*wr[gid] + Sz[gid]*ws[gid] + Tz[gid]*wt[gid];
                    
                    // Flux derivatives: ∇·F where F = u⊗u
                    // For u: ∂(u²)/∂x + ∂(uv)/∂y + ∂(uw)/∂z
                    real_t dFudx = 2*u_local*ux;
                    real_t dGudy = u_local*vy + v_local*uy;
                    real_t dHudz = u_local*wz + w_local*uz;
                    
                    // For v: ∂(uv)/∂x + ∂(v²)/∂y + ∂(vw)/∂z
                    real_t dFvdx = u_local*vx + v_local*ux;
                    real_t dGvdy = 2*v_local*vy;
                    real_t dHvdz = v_local*wz + w_local*vz;
                    
                    // For w: ∂(uw)/∂x + ∂(vw)/∂y + ∂(w²)/∂z
                    real_t dFwdx = u_local*wx + w_local*ux;
                    real_t dGwdy = v_local*wy + w_local*vy;
                    real_t dHwdz = 2*w_local*wz;
                    
                    // Volume contribution: -∇·F
                    rhsu[gid] = -(dFudx + dGudy + dHudz);
                    rhsv[gid] = -(dFvdx + dGvdy + dHvdz);
                    rhsw[gid] = -(dFwdx + dGwdy + dHwdz);
                }
            }
        }
    }
}
`

// Extract face values kernel - now for 3 components
extractFaceKernelCode := `
#define NP ` + fmt.Sprintf("%d", Np) + `
#define NFP ` + fmt.Sprintf("%d", Nfp) + `

@kernel void extractFaceValues(
    const int_t* K,
    const real_t* u_global, const int_t* u_offsets,
    const real_t* v_global, const int_t* v_offsets,
    const real_t* w_global, const int_t* w_offsets,
    real_t* Mu_global, const int_t* Mu_offsets,
    real_t* Mv_global, const int_t* Mv_offsets,
    real_t* Mw_global, const int_t* Mw_offsets,
    const int_t* VmapM
) {
    for (int part = 0; part < NPART; ++part; @outer) {
        const real_t* u = u_PART(part);
        const real_t* v = v_PART(part);
        const real_t* w = w_PART(part);
        real_t* Mu = Mu_PART(part);
        real_t* Mv = Mv_PART(part);
        real_t* Mw = Mw_PART(part);
        int k_part = K[part];
        
        // Extract values at face nodes
        for (int elem = 0; elem < KpartMax; ++elem; @inner) {
            if (elem < k_part) {
                for (int f = 0; f < 4; ++f) {
                    for (int i = 0; i < NFP; ++i) {
                        int fid = elem * NFP * 4 + f * NFP + i;
                        int vid = VmapM[f * NFP + i] + elem * NP;
                        Mu[fid] = u[vid];
                        Mv[fid] = v[vid];
                        Mw[fid] = w[vid];
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
        // 1. Extract face values to M buffers
        err = kp.RunKernelWithExtra("extractFaceValues", 
            []interface{}{vmapMDevice}, 
            "u", "v", "w", "Mu", "Mv", "Mw")
        
        // 2. Exchange partition boundaries (for all components)
        exchangePartitionBoundaries(kp, faceBuffers)
        
        // 3. Compute volume RHS (flux divergence)
        rhsuName := fmt.Sprintf("rhsu%d", stage)
        rhsvName := fmt.Sprintf("rhsv%d", stage)
        rhswName := fmt.Sprintf("rhsw%d", stage)
        err = kp.RunKernel("volumeRHS3D", 
            "u", "v", "w", 
            "ur", "us", "ut", "vr", "vs", "vt", "wr", "ws", "wt",
            "Rx", "Ry", "Rz", "Sx", "Sy", "Sz", "Tx", "Ty", "Tz", 
            rhsuName, rhsvName, rhswName)
        
        // 4. Compute flux corrections and apply LIFT
        err = kp.RunKernel("fluxCorrection3D",
            "Mu", "Mv", "Mw", "nx", "ny", "nz", "Fscale", "P_indices", "BC_types",
            "flux_u", "flux_v", "flux_w", rhsuName, rhsvName, rhswName)
        
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

- All kernels follow OCCA's @outer/@inner pattern
- Matrix operations use generated macros containing @inner loops
- MATMUL_* macros take only (IN, OUT, K_VAL) - matrix dimensions are embedded
- MATMUL_ADD_* variants accumulate to output (needed for flux corrections)
- Flux correction array allocated as device memory to avoid stack issues

### 3. Partition Parallelism

- Each partition processes independently in @outer loop
- Face exchanges occur between RHS evaluations
- Variable partition sizes handled through K array

### 4. Numerical Methods

- Conservation form with flux tensor formulation (F = **u** ⊗ **u**)
- Lax-Friedrichs numerical flux for stability
- Wall BC: reflection of normal velocity component to enforce **u**·**n̂** = 0
- LIFT matrix maps face corrections to volume contributions

### 5. Performance Considerations

- Static matrices embedded in kernel code
- Memory alignment for cache efficiency
- Coalesced access patterns in kernels
- CUDA thread limit checked during initialization

## Testing Strategy

Following the Unit Testing Principles:

1. **Start with fundamentals**: Test single element, then single partition
2. **Build systematically**: Add partitions incrementally
3. **Specific properties**: Verify conservation, convergence rates
4. **Detailed coverage**: Test boundary conditions, parallel exchange

### Test Cases

1. **Single Element Tests**
   - Constant velocity field preservation (if **u** = const, then ∂**u**/∂t = 0)
   - Linear advection accuracy
   - Divergence-free fields preservation

2. **Multi-Partition Tests**
   - Partition boundary continuity
   - Load balancing verification
   - Parallel efficiency

3. **Convergence Tests**
   - Smooth solution: O(N+1) convergence
   - Discontinuous solution: shock capturing

4. **Benchmark Problems**
   - 3D Taylor-Green vortex decay
   - Collision of vortex rings
   - Interaction of multiple vortices

## Error Handling

- Validate mesh connectivity before processing
- Check CUDA limits for partition sizes
- Verify geometric factor positivity (J > 0)
- Monitor CFL condition during time stepping
- Handle array allocation failures gracefully

## Extensions to Navier-Stokes

This Burger's equation solver demonstrates key capabilities needed for Navier-Stokes:

1. **Flux Formulation**: The conservation form with flux derivatives directly extends to the convective terms in Navier-Stokes
2. **Face Corrections**: The flux correction approach handles discontinuous solutions needed for compressible flows
3. **Partition Parallel**: The infrastructure supports the computational intensity of 3D Navier-Stokes
4. **Geometric Flexibility**: Unstructured tetrahedral meshes handle complex geometries

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

1. **Wall BC**: For the vector Burger's equation, the wall boundary condition enforces zero normal flux: **F**·**n̂** = 0. Since **F** = **u** ⊗ **u**, we have **F**·**n̂** = **u**(**u**·**n̂**). The flux correction at the wall is simply -**F**M·**n̂**, which removes the normal component of the flux.

2. **Farfield BC**: Uses characteristic analysis based on the normal velocity:
   - For outflow (**u**·**n̂** > 0): Information travels out, so we extrapolate from interior
   - For inflow (**u**·**n̂** < 0): Information comes from outside, so we impose farfield values

3. **Extension to Navier-Stokes**:
   - Wall BC will additionally need to enforce no-slip (**u** = 0) for viscous terms
   - Farfield BC will use Riemann invariants for the full system of equations (density, momentum, energy)
   - The flux tensor structure (F = **u** ⊗ **u**) extends naturally to the convective terms in Navier-Stokes

## References

- Hesthaven & Warburton, "Nodal Discontinuous Galerkin Methods"
- Spiteri & Ruuth, "Optimal Strong-Stability-Preserving Runge-Kutta Methods"
- GOCCA documentation for kernel programming patterns
- KernelProgram User's Guide for array management
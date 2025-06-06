# Approach to solving CFD/MHD in differential form with mixed elements

## Mixed Element Differential Form Implementation Summary

### Hexahedral Elements: DG-SBP with GLL Nodes

**Time advancement equation:**

```
du_i/dt = -D_ξ F_ξ - D_η F_η - D_ζ F_ζ + Σ_faces M^(-1) B_f (F* - F)_f
```

**Key components:**

- **Nodes**: Gauss-Lobatto-Legendre (GLL) points in each direction
- **Operators**: Diagonal-norm SBP operators satisfying Q = MD + D^T M = B
- **Interface flux**: SAT (Simultaneous Approximation Term) penalties
- **Norm matrix**: M = diag(w_i), where w_i are GLL quadrature weights

**Implementation details:**

- Sum-factorization for all derivative operations
- Pre-computed 1D operators, applied dimension-by-dimension
- Boundary extraction matrix B picks out face values
- Penalty parameter τ chosen for stability (typically τ ≥ 1/2)

**References:**

- Fernández et al. (2014) "Review of SBP operators for spectral methods"
- HORSES3D implementation (Ferrer et al., 2023)

### Tetrahedral Elements: Flux Reconstruction

**Time advancement equation:**

```
du_i/dt = -∇·F^D - Σ_faces ∇·[(F* - F^D) g_f]
```

Where the corrected flux is:

```
F^δ = F^D + Σ_faces (F* - F^D) g_f
```

**Key components:**

- **Nodes**: Hesthaven-Warburton warp-blend nodes for optimal conditioning
- **Polynomial basis**: Lagrange polynomials through warp-blend nodes
- **Correction function**: g_f with support on face f, zero on other faces
- **Choice of g**: Typically g_DG (recovers DG) or g_HU (energy stable)

**Implementation details:**

- Dense differentiation matrices (no tensor structure)
- Correction functions pre-computed at quadrature points
- Direct interpolation to face nodes for interface values
- Volume terms: ∇·F^D computed via differentiation matrices

**References:**

- Hesthaven & Warburton (2008) "Nodal DG Methods"
- Vincent et al. (2011) "A new class of high-order energy stable FR schemes"
- Witherden & Vincent (2014) "On the identification of symmetric quadrature
  rules for finite element methods"

### Prismatic Elements: Hybrid FR × SBP

**Time advancement equation:**

```
du_i/dt = -(D_ξ^tri ⊗ I^line)(F_x ∂ξ/∂x + F_y ∂ξ/∂y) 
          -(D_η^tri ⊗ I^line)(F_x ∂η/∂x + F_y ∂η/∂y)
          -(I^tri ⊗ D_ζ^line)F_z
          + Σ_tri_faces ∇·[(F* - F^D) g_f]      # triangular faces
          + Σ_quad_faces M^(-1) B_f (F* - F)_f   # quadrilateral faces
```

**Key components:**

- **Nodes**: Hesthaven-Warburton (triangle) × GLL (line direction)
- **Triangle operators**: FR differentiation matrices D_ξ^tri, D_η^tri
- **Line operator**: 1D SBP operator D_ζ^line with diagonal norm
- **Mixed face treatment**: FR corrections on triangular faces, SAT on quad
  faces

**Implementation details:**

- Partial tensor-product structure enables semi-sum-factorization
- Triangular caps use scalar correction function g
- Quadrilateral sides use tensor-product SAT
- Metric terms (∂ξ/∂x, etc.) computed via reference element mapping

**References:**

- Williams & Jameson (2013) "Energy stable FR schemes for mixed elements"
- NodesAndModes.jl for prism node generation

### Pyramid Elements: Collapsed Coordinate FR

**Time advancement equation:**

```
du_i/dt = -∇·F^D - ∇·[(F*_tri - F^D_tri) g_tri] 
          - Σ_quad_faces ∇·[(F*_quad - F^D_quad) g_quad]
```

**Key components:**

- **Nodes**: Collapsed coordinate mapping from hex reference element
- **Special treatment**: Singularity at apex requires careful differentiation
- **Mixed faces**: One quad base, four triangular sides

**Note**: Pyramids are the most challenging element type with limited high-order
implementations in production codes.

### Unified Mortar Coupling

All element types couple through the mortar interface via:

1. **Project to mortar**: u_mortar = Π_h(u_elem) via L² projection
2. **Compute flux**: F*_mortar = Riemann(u_L, u_R, n)
3. **Project back**: Each element receives F* and applies its correction:
    - Hex: SAT penalty
    - Tet/Tri faces: FR correction g
    - Quad faces: Tensor SAT

The key insight is that the mortar doesn't need to know about element-specific
correction mechanisms - it just provides the common flux F*.

**Conservation**: Maintained through the discrete compatibility condition:

```
∫_Γ F*_mortar · v_A dΓ = -∫_Γ F*_mortar · v_B dΓ
```

This holds regardless of element type or correction method used.

## Question: Is there one unified time advancement?

## Unified time advancement structure:

```python
# Single time advancement loop for ALL elements
for stage in RK_stages:
    # Compute residuals for all elements
    for elem in all_elements:
        if elem.type == HEX:
            residual[elem] = compute_hex_residual(u[elem])
        elif elem.type == TET:
            residual[elem] = compute_tet_residual(u[elem])
        elif elem.type == PRISM:
            residual[elem] = compute_prism_residual(u[elem])
    
    # Single unified update
    u_new = u + dt * residual
```

## Why element-specific residual functions?

Each element type computes its spatial derivatives differently:

- **Hex**: Sum-factorization with SBP operators
- **Tet**: Dense differentiation matrices with FR correction
- **Prism**: Hybrid tensor/non-tensor operations

But they all produce the same output: `-∇·F + source` evaluated at the element's
nodes.

## The time advancement is identical:

Whether using RK4, SSPRK3, or any other time integrator, you advance ALL
elements together:

```
k1 = f(u^n, t^n)
k2 = f(u^n + dt/2 * k1, t^n + dt/2)
...
u^(n+1) = u^n + dt/6 * (k1 + 2k2 + 2k3 + k4)
```

The function `f()` calls the appropriate residual computation for each element
type, but the time integration treats all DOFs uniformly.

## Cache efficiency consideration:

You mentioned grouping like elements - this happens at the **residual
computation level**:

```python
# Cache-efficient residual computation
def compute_all_residuals(u):
    # Process all hexes together
    for elem in hex_elements:
        residual[elem] = compute_hex_residual(u[elem])
    
    # Process all tets together  
    for elem in tet_elements:
        residual[elem] = compute_tet_residual(u[elem])
    
    # Process all prisms together
    for elem in prism_elements:
        residual[elem] = compute_prism_residual(u[elem])
    
    return residual
```

So: **one time advancement equation** (du/dt = residual), but **element-specific
residual implementations** that are grouped by type for efficiency.

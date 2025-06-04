# Discontinuous Galerkin Face Quadrature and Lift Operator Implementation Guide

## Overview

This document provides detailed implementation instructions for face quadrature points and lift operator construction in a 3D discontinuous Galerkin element library. The lift operator maps face flux contributions to volume residuals, ensuring conservative flux transfer between elements.

## Mathematical Foundation

The lift operator L transforms surface integrals to volume contributions via the relation:
```
∫_K L[f] · φ_i dx = ∫_∂K f · φ_i ds
```

In matrix form: **L = M⁻¹B^T M_f**, where:
- M is the element mass matrix
- B is the face-to-volume interpolation operator
- M_f is the face mass matrix

For conservation, the numerical flux must be unique at interfaces, and the lift operator must accurately integrate this flux using consistent quadrature rules on both sides of each face.

**Key References:**
- Hesthaven & Warburton (2008), "Nodal Discontinuous Galerkin Methods", Chapters 6-7
- Kopriva (2009), "Implementing Spectral Methods for PDEs", Section 7.3
- Cockburn & Shu (1998), "The Local Discontinuous Galerkin Method", Section 3

---

## Tetrahedral Elements (PKD Basis)

### Step 1: Face Quadrature Point Definition

Each tetrahedral element has 4 triangular faces. For polynomial order p, use symmetric quadrature rules on triangles.

**Implementation Steps:**
1. Define reference triangle coordinates: (0,0), (1,0), (0,1)
2. Load symmetric quadrature rules for triangles (Dunavant 1985 or Wandzura-Xiao 2003)
3. Store quadrature points ξ_q, η_q and weights w_q for order 2p+1 accuracy

**Unit Tests:**
- Verify quadrature integrates monomials r^i s^j exactly for i+j ≤ 2p+1
- Check sum of weights equals triangle area (0.5)
- Validate all quadrature points lie within reference triangle

### Step 2: Face-to-Volume Mapping

Map 2D face quadrature points to 3D tetrahedral reference coordinates.

**Implementation Steps:**
1. Define face-to-tet mappings for each of 4 faces:
    - Face 0 (vertices 1,2,3): (ξ,η) → (ξ, η, 1-ξ-η)
    - Face 1 (vertices 0,2,3): (ξ,η) → (0, η, 1-η)
    - Face 2 (vertices 0,1,3): (ξ,η) → (ξ, 0, 1-ξ)
    - Face 3 (vertices 0,1,2): (ξ,η) → (ξ, η, 0)
2. Compute outward normal vectors for each face in reference coordinates

**Unit Tests:**
- Verify mapped points lie on correct face
- Check normal vectors have unit length and correct orientation
- Validate face vertices map correctly

### Step 3: Basis Function Evaluation at Face Quadrature Points

Evaluate PKD basis functions at all face quadrature points.

**Implementation Steps:**
1. For each face f and quadrature point q:
    - Map (ξ_q, η_q) to 3D reference coordinates
    - Evaluate all PKD basis functions φ_i at this point
    - Store in matrix B[f,i,q] = φ_i(x_q)
2. Include Jacobian scaling for physical elements

**Unit Tests:**
- Verify basis functions sum to 1 at each quadrature point
- Check polynomial reproduction on faces
- Validate against analytical PKD formulas at vertices

### Step 4: Face Mass Matrix Construction

Build face mass matrices for flux integration.

**Implementation Steps:**
1. For each face f:
    - M_f[i,j] = Σ_q w_q φ_i(x_q) φ_j(x_q) |J_f|
    - J_f is the face Jacobian (area scaling)
2. Store as dense matrices (small size for faces)

**Unit Tests:**
- Verify symmetry and positive definiteness
- Check conditioning number is reasonable
- Validate mass conservation: row sums equal face area

### Step 5: Lift Operator Assembly

Construct the complete lift operator for each face.

**Implementation Steps:**
1. Compute element mass matrix M (if not already available)
2. For each face f:
    - Form B_f^T M_f where B_f contains basis evaluations
    - Solve M L_f^T = B_f^T M_f for L_f
    - Store L_f as dense matrix or in compressed format
3. Optionally combine into single lift operator

**Unit Tests:**
- Apply lift to constant flux, verify volume integral equals surface integral
- Check polynomial flux lifting preserves accuracy
- Validate conservation: sum of lifted values equals face integral

---

## Hexahedral Elements (Tensor Product Basis)

### Step 1: Face Quadrature Point Definition

Each hexahedral element has 6 quadrilateral faces. Use tensor product Gauss-Legendre quadrature.

**Implementation Steps:**
1. Generate 1D Gauss-Legendre points and weights for order p+1
2. Form 2D tensor product quadrature on [-1,1]²
3. Store (ξ_i, η_j) points and w_i × w_j weights

**Unit Tests:**
- Verify exact integration of x^i y^j for i,j ≤ 2p+1
- Check weight sum equals 4 (area of reference square)
- Validate quadrature point symmetry

### Step 2: Face-to-Volume Mapping

Map 2D face quadrature to 3D hex reference coordinates [-1,1]³.

**Implementation Steps:**
1. Define mappings for 6 faces:
    - Face 0 (ξ=-1): (η,ζ) → (-1,η,ζ)
    - Face 1 (ξ=+1): (η,ζ) → (+1,η,ζ)
    - Face 2 (η=-1): (ξ,ζ) → (ξ,-1,ζ)
    - Face 3 (η=+1): (ξ,ζ) → (ξ,+1,ζ)
    - Face 4 (ζ=-1): (ξ,η) → (ξ,η,-1)
    - Face 5 (ζ=+1): (ξ,η) → (ξ,η,+1)
2. Store constant outward normals for each face

**Unit Tests:**
- Verify face orientation consistency
- Check point mapping preserves tensor structure
- Validate normal vectors sum to zero

### Step 3: Tensor Product Basis Evaluation

Efficiently evaluate tensor product basis at face points.

**Implementation Steps:**
1. Use 1D basis evaluations to build face evaluations:
    - For ξ-faces: φ(x_q) = φ_ξ(±1) × φ_η(η_q) × φ_ζ(ζ_q)
    - Similar for η and ζ faces
2. Exploit tensor structure for efficiency
3. Store sparse patterns (many zeros due to tensor structure)

**Unit Tests:**
- Check basis partition of unity on faces
- Verify tensor product structure preserved
- Validate against direct 3D evaluation

### Step 4: Structured Face Mass Matrix

Exploit tensor product structure for face mass matrices.

**Implementation Steps:**
1. Face mass matrices are tensor products of 1D mass and identity:
    - ξ-faces: M_f = I ⊗ M_1D ⊗ M_1D
    - Similar patterns for other face orientations
2. Store only 1D components for efficiency

**Unit Tests:**
- Verify Kronecker product structure
- Check eigenvalue patterns match theory
- Validate integration accuracy

### Step 5: Efficient Lift Operation

Use sum factorization for efficient lift application.

**Implementation Steps:**
1. Implement sum-factorization lift without forming full matrices
2. For ξ-faces: apply 1D operations in sequence:
    - Contract in η direction
    - Contract in ζ direction
    - Scale by inverse mass in ξ
3. Store only 1D operators

**Unit Tests:**
- Compare against dense lift operation
- Verify O(p^4) vs O(p^6) scaling
- Check conservation properties maintained

---

## Prismatic Elements (Hybrid Basis)

### Step 1: Mixed Face Types

Prisms have 3 quadrilateral faces and 2 triangular faces.

**Implementation Steps:**
1. For quad faces: use tensor product Gauss-Legendre quadrature
2. For tri faces: use symmetric triangle quadrature (same as tets)
3. Store separate quadrature data for each face type

**Unit Tests:**
- Verify appropriate accuracy for each face type
- Check quadrature weight consistency
- Validate point distributions

### Step 2: Collapsed Coordinate Mapping

Handle the triangular-to-quadrilateral collapse carefully.

**Implementation Steps:**
1. Quad face mappings (sides):
    - Faces perpendicular to triangular cross-section
    - Standard tensor product mappings
2. Tri face mappings (top/bottom):
    - Use collapsed coordinate transformation
    - Account for coordinate singularity at collapse edge
3. Compute metric terms accounting for collapse

**Unit Tests:**
- Check Jacobian remains positive
- Verify smooth variation near collapsed edge
- Validate normal vector computation

### Step 3: Hybrid Basis Evaluation

Evaluate tensor×simplex basis at face quadrature points.

**Implementation Steps:**
1. For quad faces: tensor product of 1D × triangle basis
2. For tri faces: pure 2D PKD basis evaluation
3. Handle collapsed coordinate carefully for stability

**Unit Tests:**
- Verify basis completeness on each face type
- Check continuity across edges
- Validate polynomial accuracy preservation

### Step 4: Face-Specific Mass Matrices

Build mass matrices respecting face geometry.

**Implementation Steps:**
1. Quad faces: exploit partial tensor structure
2. Tri faces: standard symmetric mass matrix
3. Include appropriate metric scaling

**Unit Tests:**
- Check mass conservation per face type
- Verify expected sparsity patterns
- Validate conditioning

### Step 5: Unified Lift Strategy

Combine different face types into single lift framework.

**Implementation Steps:**
1. Maintain separate lift operators per face
2. Apply appropriate lift based on face type
3. Ensure consistent flux integration accuracy

**Unit Tests:**
- Test mixed face type configurations
- Verify conservation across face types
- Check accuracy for smooth solutions

---

## Pyramid Elements

### Step 1: Complex Face Configuration

Pyramids have 1 quadrilateral base and 4 triangular faces.

**Implementation Steps:**
1. Quad base: tensor product quadrature
2. Triangular faces: symmetric quadrature
3. Account for varying face sizes in physical space

**Unit Tests:**
- Verify quadrature accuracy per face
- Check geometric consistency
- Validate reference mappings

### Step 2: Rational Basis Handling

Pyramid bases often involve rational functions.

**Implementation Steps:**
1. Carefully evaluate rational basis at face points
2. Handle potential singularities at apex
3. Use stable evaluation algorithms (Bergot & Duruflé 2010)

**Unit Tests:**
- Check numerical stability near apex
- Verify polynomial reproduction
- Test conditioning of evaluations

### Step 3: Apex Singularity Management

Special treatment needed near pyramid apex.

**Implementation Steps:**
1. Use regularized coordinates near apex
2. Implement stable quadrature for triangular faces
3. Ensure continuous flux approximation

**Unit Tests:**
- Test flux conservation near apex
- Verify stability for high-order approximations
- Check error convergence rates

### Step 4: Mixed Mass Matrix Assembly

Different strategies for base vs side faces.

**Implementation Steps:**
1. Base: exploit tensor product structure
2. Sides: account for varying triangle shapes
3. Include geometric factors correctly

**Unit Tests:**
- Verify mass conservation
- Check symmetry properties
- Test scaling with element distortion

### Step 5: Pyramid-Specific Lift Operator

Account for unique pyramid geometry in lift construction.

**Implementation Steps:**
1. Separate treatment for base and side faces
2. Careful handling of apex-adjacent faces
3. Ensure conservation despite rational basis

**Unit Tests:**
- Verify lift operator accuracy
- Test conservation properties
- Check stability for distorted pyramids

---

## Universal Mortar Implementation

### Overview

The mortar provides a common space for flux computation between potentially non-conforming elements of different types.

### Step 1: Mortar Space Definition

Define polynomial space on each face type that accommodates all adjacent elements.

**Implementation Steps:**
1. For each face geometry (tri/quad):
    - Define maximum polynomial order p_max
    - Choose quadrature order ≥ 2p_max + 1
2. Use consistent orientation convention
3. Store mortar basis functions

**Unit Tests:**
- Verify polynomial completeness
- Check quadrature accuracy
- Validate basis orthogonality

### Step 2: Element-to-Mortar Projection

Build projection operators from each element type to mortar space.

**Implementation Steps:**
1. For each element/face combination:
    - Evaluate element basis at mortar quadrature points
    - Form L² projection matrix
    - P_mortar = M_mortar⁻¹ B^T W
2. Account for different element orientations
3. Precompute and store projections

**Unit Tests:**
- Verify projection preserves constants
- Check accuracy for polynomials up to p
- Test orthogonality of projection

### Step 3: Flux Computation Infrastructure

Enable efficient flux evaluation on mortar space.

**Implementation Steps:**
1. Structure for storing left/right states on mortar
2. Interface for Riemann solver evaluation
3. Vectorized flux computation at all quadrature points

**Unit Tests:**
- Test simple flux functions (upwind)
- Verify vectorization correctness
- Check flux consistency

### Step 4: Mortar-to-Element Distribution

Project fluxes back from mortar to element faces.

**Implementation Steps:**
1. Build distribution operators (transpose of projection)
2. Ensure conservation: Σ(element fluxes) = 0
3. Include Jacobian scaling for physical elements

**Unit Tests:**
- Verify conservation to machine precision
- Test flux distribution accuracy
- Check scaling with element size

### Step 5: Non-Conforming Interface Handling

Extend mortar for h/p non-conforming interfaces.

**Implementation Steps:**
1. Detect non-conforming interfaces
2. Build appropriate sub-mortars for h-refinement
3. Use higher-order mortar for p-refinement
4. Maintain strict conservation

**Unit Tests:**
- Test 2:1 h-refinement conservation
- Verify p-refinement accuracy
- Check complex interface configurations

---

## Implementation Notes

1. **Memory Layout**: Store face data contiguously for cache efficiency
2. **Precomputation**: All geometric quantities computed during initialization
3. **Parallelization**: Face computations naturally parallel
4. **Precision**: Use appropriate quadrature precision for conservation

## References

- Bergot, M. & Duruflé, M. (2010). "High-order optimal edge elements for pyramids, prisms and hexahedra". J. Comput. Phys.
- Cockburn, B. & Shu, C.W. (1998). "The Local Discontinuous Galerkin Method for Time-Dependent Convection-Diffusion Systems". SIAM J. Numer. Anal.
- Dunavant, D.A. (1985). "High degree efficient symmetrical Gaussian quadrature rules for the triangle". Int. J. Numer. Methods Eng.
- Hesthaven, J.S. & Warburton, T. (2008). "Nodal Discontinuous Galerkin Methods". Springer.
- Kopriva, D.A. (2009). "Implementing Spectral Methods for Partial Differential Equations". Springer.
- Wandzura, S. & Xiao, H. (2003). "Symmetric quadrature rules on a triangle". Comput. Math. Appl.
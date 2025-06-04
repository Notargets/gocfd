# Symmetric Quadrature Rules for 3D Finite Elements: A Comprehensive Technical Report

## Advanced quadrature formulations for tetrahedral, hexahedral, prismatic, and pyramidal elements

This report presents a comprehensive analysis of symmetric quadrature rules for 3D finite elements, focusing on Williams-Shunn-Jameson rules for tetrahedra and optimal collocated rules for other element types that enable exact integration of polynomials up to order 2P for order P elements.

## Williams-Shunn-Jameson Rules for Tetrahedra

The Williams-Shunn-Jameson symmetric quadrature rules represent a significant advancement in tetrahedral numerical integration, building upon sphere close-packed (SCP) lattice arrangements to achieve optimal symmetric point distributions. These rules extend the earlier work of Shunn and Ham (2012) through advanced optimization techniques.

### Mathematical Foundation and Properties

The core innovation lies in using cubic close-packed lattice arrangements that preserve tetrahedral symmetry under affine transformations. The optimization process determines point configurations and weights that maintain full tetrahedral symmetry while minimizing the number of quadrature points for a given polynomial degree.

**Available Rules and Specifications:**
- **Order 1**: 1 point (centroid), exact for degree 1
- **Order 2-3**: 4 points (vertices), exact for degrees 2-3
- **Order 4-5**: 10 points (CCP lattice), exact up to degree 5
- **Order 6**: 20 points (extended CCP), exact up to degree 6
- **Order 7**: 35 points (higher-order CCP), exact up to degree 7
- **Order 8**: 56 points (advanced CCP), exact up to degree 8
- **Order 9**: 84 points (Williams-Shunn-Jameson extension), exact up to degree 9

The rules achieve **subgeometric exponential convergence** with index values 0.5 < r < 1, and all weights remain positive. Points are organized by symmetry groups: vertex points, edge points, face points, and interior points.

### Key Citations and Implementation

**Primary Sources:**
1. Williams, D.M., Shunn, L., Jameson, A. (2014). "Symmetric quadrature rules for simplexes based on sphere close packed lattice arrangements." *Journal of Computational and Applied Mathematics*, 266, 18–38.
2. Shunn, L., Ham, F. (2012). "Symmetric quadrature rules for tetrahedra based on a cubic close-packed lattice arrangement." *Journal of Computational and Applied Mathematics*, 236(17), 4348–4364.

**Implementation:** The QuadPy Python library provides the most accessible implementation through `quadpy.t3._williams_shunn_jameson` modules, with functions like `quadpy.t3.williams_shunn_jameson_*()` for different orders.

## Hexahedral Elements with Colocation Property

For hexahedral elements, **Gauss-Lobatto-Legendre (GLL) quadrature** provides the optimal solution for achieving colocation with nodal points while maintaining exact integration for polynomials up to order 2P-1 for order P elements.

### GLL Quadrature Specifications

For polynomial order P, the GLL quadrature uses N = P+1 points and achieves exact integration for polynomials up to degree 2N-3 = 2P-1. The critical advantage is that quadrature points include the element boundaries (-1, +1), enabling natural boundary condition enforcement and element connectivity.

**1D GLL Points and Weights:**

**P=1 (N=2):**
- Points: x₁ = -1, x₂ = +1
- Weights: w₁ = 1, w₂ = 1

**P=2 (N=3):**
- Points: x₁ = -1, x₂ = 0, x₃ = +1
- Weights: w₁ = 1/3, w₂ = 4/3, w₃ = 1/3

**P=3 (N=4):**
- Points: x₁ = -1, x₂ = -√(1/5), x₃ = √(1/5), x₄ = +1
- Weights: w₁ = 1/6, w₂ = 5/6, w₃ = 5/6, w₄ = 1/6

**P=4 (N=5):**
- Points: x₁ = -1, x₂ = -√(3/7), x₃ = 0, x₄ = √(3/7), x₅ = +1
- Weights: w₁ = 1/10, w₂ = 49/90, w₃ = 32/45, w₄ = 49/90, w₅ = 1/10

For 3D hexahedral elements, the quadrature uses tensor products:
- Total points: N³
- Coordinates: (ξᵢ, ηⱼ, ζₖ) where each is a 1D GLL point
- Weights: wᵢⱼₖ = wᵢ × wⱼ × wₖ

### Mathematical Properties and Citations

The colocation property creates diagonal mass matrices through the summation-by-parts property, enabling efficient explicit time-stepping. This property is mathematically proven through:
∑ᵢ wᵢ ∂φⱼ/∂ξ|ξᵢ φₖ(ξᵢ) = δⱼₖ

**Key References:**
1. Kopriva, D.A., Gassner, G. (2010). "On the Quadrature and Weak Form Choices in Collocation Type Discontinuous Galerkin Spectral Element Methods." *J. Sci. Comput.* 44, 136–155.
2. Gassner, G., Kopriva, D.A. (2011). "A comparison of the dispersion and dissipation errors of Gauss and Gauss–Lobatto discontinuous Galerkin spectral element methods." *SIAM J. Sci. Comput.* 33, 2560–2579.

## Prismatic Elements: Efficient Non-Product Rules

Prismatic elements benefit from specialized non-product quadrature rules that achieve approximately 50% reduction in quadrature points compared to tensor-product approaches.

### Optimal Quadrature Specifications

**Reference Domain:** Unit triangular prism with triangle base {(ξ,η) | ξ,η ≥ 0, ξ+η ≤ 1} and height ζ ∈ [-1,1].

**Efficient Non-Product Rules (Wandzura & Xiao):**
- **Degree 2**: 6 points (vs 9 for tensor product)
- **Degree 3**: 8 points (vs 12 for tensor product)
- **Degree 5**: 21 points (vs 35 for tensor product)
- **Degree 9**: 73 points (achieves theoretical minimum)

These rules use polynomial moment fitting with symmetric point arrangements to minimize the number of quadrature points while maintaining accuracy.

**Key Citation:**
Wandzura, S., Xiao, H. (2003). "New computationally efficient quadrature formulas for triangular prism elements." *Engineering Computations* 20(3), 390-412.

## Pyramidal Elements: Handling the Apex Singularity

Pyramidal elements present unique challenges due to the geometric singularity at the apex where the quadrilateral base contracts to a point. The breakthrough solution uses Jacobi polynomials with appropriate weight functions.

### Bergot-Cohen-Duruflé Method

The method employs Jacobi polynomials orthogonal with respect to weight (1-z)² to handle the singularity. The basis functions use polynomials up to order r in x,y coordinates but non-polynomial rational functions in the z-coordinate.

**Reference Pyramid:** {(x,y,z) | |x|,|y| ≤ 1-z, 0 ≤ z ≤ 1}

**Low-Order Quadrature Data:**

**Order 1 (1 point):**
- Point: (0, 0, 0.25)
- Weight: 8/6
- Polynomial degree: 1

**Order 2 (5 points):**
- Vertex configuration with positive weights
- Exact for polynomials up to degree 2

**Implementation Approach:**
1. Use conical product rules: Legendre × Legendre × Jacobi(2,0)
2. The Jacobi rule automatically includes the (1-z)² factor
3. Points arranged on parallel planes z = constant for stability

**Key Citations:**
1. Bergot, M., Cohen, G., Duruflé, M. (2010). "Higher-order finite elements for hybrid meshes using new nodal pyramidal elements." *J. Sci. Comput.* 42, 345-381.
2. Witherden, F.D., Vincent, P.E. (2015). "On the identification of symmetric quadrature rules for finite element methods." *Computers & Mathematics with Applications* 69, 1232-1241.

## Software Implementation and Availability

### Primary Software Resources

**QuadPy (Python):** Comprehensive library containing Williams-Shunn-Jameson rules for tetrahedra and various rules for other elements. Access via `pip install quadpy`.

**Polyquad (C++):** Witherden-Vincent implementation for systematic generation of symmetric quadrature rules with positive weights for all element types.

**HORSES3D:** Open-source high-order DG solver implementing GLL quadrature for hexahedral elements with arbitrary polynomial orders.

**Nektar++:** Spectral/hp element framework supporting multiple quadrature types with extensive documentation.

### Integration Guidelines

For implementation, ensure:
1. Quadrature order ≥ 2P for order P elements
2. Validate polynomial exactness through systematic testing
3. Monitor condition numbers for high-order implementations
4. Use appropriate transformations for reference-to-physical element mapping

## Comparative Analysis and Recommendations

### Element-Specific Recommendations

**Tetrahedra:** Use Williams-Shunn-Jameson rules up to order 9; for higher orders, consider Jaśkowiec-Sukumar rules (up to degree 20).

**Hexahedra:** GLL quadrature provides optimal balance between accuracy (2P-1) and functionality (colocation, boundary points).

**Prisms:** Implement Wandzura-Xiao non-product rules for computational efficiency; use tensor products only for rapid prototyping.

**Pyramids:** Apply Bergot method with Jacobi polynomial approach; use rational basis functions to handle apex singularity properly.

### Performance Considerations

The colocation property in all element types enables:
- Diagonal mass matrices for efficient time-stepping
- Natural boundary condition enforcement
- Element-to-element connectivity
- Reduced memory requirements

Trade-offs include potentially more restrictive time-step constraints and the need for specialized basis functions (particularly for pyramids).

## Conclusion

This comprehensive analysis demonstrates that optimal symmetric quadrature rules exist for all standard 3D finite element types, each addressing specific mathematical and computational challenges. The Williams-Shunn-Jameson rules excel for tetrahedral elements through sphere close-packed lattice optimization, while GLL quadrature provides the ideal collocated solution for hexahedra. Prismatic elements benefit from efficient non-product rules, and pyramidal elements require sophisticated Jacobi polynomial approaches to handle their inherent singularity. Together, these quadrature formulations enable robust, high-order finite element methods on general polyhedral meshes.
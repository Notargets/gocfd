# Williams-Shunn-Jameson symmetric quadrature rules for tetrahedral integration

## Core findings and polynomial integration capability

The Williams-Shunn-Jameson symmetric quadrature rules provide highly efficient numerical integration on tetrahedra, with each order exactly integrating polynomials up to its corresponding degree. **Order 1 integrates degree 1 polynomials with 1 point, order 2 handles degree 2 with 4 points, order 3 manages degree 3 with 10 points, order 4 covers degree 4 with 20 points, and order 5 integrates degree 5 polynomials using 35 points**. These rules originate from two seminal papers: Shunn & Ham (2012) establishing the foundation using cubic close-packed lattice arrangements, and Williams, Shunn & Jameson (2014) extending the methodology to higher orders and multiple dimensions.

## Reference tetrahedron and coordinate systems

The quadrature rules use the **unit simplex tetrahedron** as their reference element, with vertices at **(0,0,0), (1,0,0), (0,1,0), and (0,0,1)**. Points within this tetrahedron are expressed using barycentric coordinates (ξ₁, ξ₂, ξ₃, ξ₄) where ξ₁ + ξ₂ + ξ₃ + ξ₄ = 1 and all ξᵢ ≥ 0. The unit simplex has a volume of 1/6, which is crucial for weight calculations. All quadrature points are strictly interior to the tetrahedron, and all weights are positive, ensuring numerical stability.

## Transformation to [-1,1]³ reference tetrahedron

To map these rules to a tetrahedron defined over [-1,1]³ coordinates, you need specific transformation formulas. For a common [-1,1]³ reference tetrahedron with vertices at **(-1,-1,-1), (1,-1,-1), (-1,1,-1), and (-1,-1,1)**, the transformation from unit simplex coordinates (ξ₁, ξ₂, ξ₃, ξ₄) to (r,s,t) coordinates follows:

```
r = 2ξ₁ - 1
s = 2ξ₂ - 1  
t = 2ξ₃ - 1
```

The inverse transformation from (r,s,t) to barycentric coordinates:
```
ξ₁ = (1 + r)/2
ξ₂ = (1 + s)/2
ξ₃ = (1 + t)/2
ξ₄ = 1 - ξ₁ - ξ₂ - ξ₃ = (1 - r - s - t)/2
```

**The Jacobian determinant for this transformation is 8**, meaning quadrature weights must be multiplied by 8 when transforming from the unit simplex (volume 1/6) to the [-1,1]³ reference tetrahedron (volume 4/3).

## Available numerical values for order 4 (20-point rule)

Our research uncovered exact values for the 20-point rule (order 4), which consists of three symmetry groups:

**Group 1 (4 points)** - One coordinate dominant:
- Barycentric: [0.03235259, 0.03235259, 0.90294222, 0.03235259] and cyclic permutations
- Weight: 0.00706707 per point

**Group 2 (12 points)** - Two coordinates balanced:
- Barycentric: [0.06036044, 0.26268258, 0.61659653, 0.06036044] and permutations
- Weight: 0.04699867 per point

**Group 3 (4 points)** - Three coordinates equal:
- Barycentric: [0.3097693, 0.3097693, 0.07069209, 0.3097693] and permutations
- Weight: 0.04699867 per point

To generate all 20 points, apply all unique permutations preserving the symmetry pattern of each group.

## Obtaining complete numerical values for orders 1-5

The complete high-precision numerical values for all orders are available through several sources:

**1. QuadPy Python library** (most accessible):
```python
import quadpy
for order in range(1, 6):
    scheme = quadpy.tetrahedron.ShunnHam(order)
    print(f"Order {order}: {len(scheme.points)} points")
    print("Points:", scheme.points)  # In barycentric coordinates
    print("Weights:", scheme.weights)
```

**2. Original academic papers**:
- Shunn & Ham (2012) contains complete tables for orders 1-6 in the appendices
- Williams, Shunn & Jameson (2014) extends to higher orders with 84-point rules

**3. John Burkardt's scientific computing libraries** provide FORTRAN90, C++, and MATLAB implementations with embedded numerical values, though these may include Keast or Felippa rules rather than Williams-Shunn-Jameson specifically.

## Implementation strategy for golang

For implementing these rules in golang, the recommended approach involves:

1. **Store quadrature data** in a structure mapping order to points and weights
2. **Use barycentric coordinates** internally for symmetry preservation
3. **Implement transformation functions** to convert between coordinate systems
4. **Apply the Jacobian scaling** when transforming to different reference elements
5. **Validate implementation** by integrating monomials up to the specified degree

The quadrature rules' **symmetric properties** mean you can store only the unique representatives and generate full point sets through permutation. This reduces storage requirements while maintaining numerical precision.

## Critical implementation considerations

These quadrature rules achieve remarkable efficiency through their symmetric arrangement based on sphere close-packed lattices. **All weights remain positive**, preventing sign cancellation errors, and **all points lie strictly inside the tetrahedron**, avoiding boundary singularities. The rules exhibit optimal or near-optimal efficiency for their polynomial degrees, making them particularly valuable for high-order finite element methods.

When implementing transformations, remember that the **volume scaling factor** between reference elements affects the quadrature weights. For arbitrary tetrahedra, use affine transformations with the Jacobian determinant to scale weights appropriately. The high precision of the original values (often 15-34 decimal places) should be preserved to maintain accuracy in high-order methods.
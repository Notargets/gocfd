# Hierarchical and Incremental, Specific and Detailed Unit Testing Principles

## Core Testing Philosophy

Unit tests should validate mathematical and algorithmic correctness through
systematic, progressive verification that builds confidence from fundamental
properties to complex behaviors.

## 1. Hierarchical Testing

- **Start with fundamentals**: Test the simplest cases first (constants, linear
  functions, low-order polynomials)
- **Build systematically**: Each level of complexity should build upon validated
  lower levels
- **Preserve properties**: Verify that properties proven at order N-1 remain
  valid at order N
- **Example**: For polynomial interpolation, test constant → linear →
  quadratic → cubic, ensuring each level maintains exactness of all lower levels

## 2. Incremental Validation

- **Progressive complexity**: Test across increasing orders/dimensions/sizes
  systematically
- **No gaps**: Test every intermediate level, not just extremes (test
  N=1,2,3,4,5,6, not just N=1,6)
- **Cross-validation**: Properties valid for specific cases should hold for
  general cases
- **Example**: If differentiation is exact for polynomials of degree p at order
  N, verify this remains true at orders N+1, N+2, etc.

## 3. Specific Property Testing

- **Mathematical exactness**: Test specific mathematical properties that must
  hold exactly (not approximately)
- **Known solutions**: Use problems with analytical solutions to verify
  numerical accuracy
- **Orthogonality/Identity**: Test fundamental mathematical relationships (V *
  V^(-1) = I, orthogonal basis properties)
- **Conservation laws**: Verify that conserved quantities remain constant
- **Example**: Derivative of x³ must be exactly 3x², not approximately

## 4. Detailed Coverage

- **Complete basis testing**: For polynomial methods, test ALL monomials x^i *
  y^j * z^k up to the design order
- **Boundary conditions**: Test behavior at element boundaries, interfaces, and
  special points
- **Degeneracies**: Test special cases (zero values, repeated nodes, aligned
  points)
- **Numerical stability**: Verify condition numbers and error propagation across
  orders

## 5. Implementation Guidelines

- **Test independence**: Each test should validate one specific property
- **Clear test names**: Name tests to indicate exactly what property is being
  validated
- **Quantitative thresholds**: Use appropriate tolerances (1e-10 for exact math,
  1e-8 for floating point)
- **Progressive test structure**: Order tests from simple to complex within test
  suites
- **Reusable test infrastructure**: Create helper functions for common
  validation patterns

## 6. Error Analysis

- **Expected convergence rates**: Verify that errors decrease at the theoretical
  rate with increasing order
- **Error accumulation**: Test how errors propagate through multiple operations
- **Relative vs absolute error**: Use appropriate error metrics for the quantity
  being tested

## Example Test Progression

For a spectral element method of order N:

1. Test nodes are correctly positioned (no duplicates, within domain)
2. Test Vandermonde matrix is invertible
3. Test constant function interpolation (exact)
4. Test linear function interpolation (exact)
5. Test polynomial interpolation up to degree N (exact)
6. Test differentiation of constants (exactly zero)
7. Test differentiation of linear functions (exact)
8. Test differentiation of polynomials up to degree N-1 (exact)
9. Test that all above properties hold for N=1 through N=6

## Key Principle

**"If it works for order N-1, it must still work for order N"** - Higher-order
methods should maintain all the accuracy properties of lower-order methods while
adding new capabilities.
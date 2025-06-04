# 3D Finite Element Polynomial Basis Implementation Plan

## Phase 1: Tetrahedral Elements (PKD Basis)

### Step 1.1: Koornwinder Basis Functions with Dubiner Normalization
**Implementation:**
```go
// koornwinder.go
type TetBasis struct {
    Order int
    Np    int // Number of basis functions = (P+1)(P+2)(P+3)/6
}

// Koornwinder-Dubiner orthonormal basis on reference tetrahedron
// Reference: Hesthaven & Warburton (2008), "Nodal Discontinuous Galerkin Methods", Ch. 10
func (tb *TetBasis) EvaluateBasis(r, s, t float64) []float64
func (tb *TetBasis) ComputeVandermonde(r, s, t []float64) [][]float64

// Collapsed coordinate transform
func CollapseCoordinates(r, s, t float64) (a, b, c float64) {
    // a = 2(1+r)/(1-s-t) - 1 if s+t != 1, else a = -1
    // b = 2(1+s)/(1-t) - 1 if t != 1, else b = -1  
    // c = t
}

// Normalized Jacobi polynomials P_n^(alpha,beta)(x)
func JacobiP(n int, alpha, beta, x float64) float64

// PKD basis formula: ψ_ijk(r,s,t) = P_i^(0,0)(a) * ((1-b)/2)^i * P_j^(2i+1,0)(b) * ((1-c)/2)^(i+j) * P_k^(2i+2j+2,0)(c) * N_ijk
// where N_ijk is the normalization factor for orthonormality
func (tb *TetBasis) ComputeNormalizationFactor(i, j, k int) float64
```

**Unit Tests:**
```go
// Test 1.1.1: Basis orthonormality verification
func TestKoornwinderOrthonormality(t *testing.T)
// - Integrate ψ_p * ψ_q over reference tetrahedron
// - Verify <ψ_p, ψ_q> = δ_pq within tolerance 1e-12
// - Use high-order Grundmann-Moeller cubature

// Test 1.1.2: Polynomial reproduction
func TestPolynomialInterpolation(t *testing.T)
// - Interpolate monomials r^i * s^j * t^k for i+j+k <= P
// - Verify exact reproduction at Warburton nodes
// - Check condition number of Vandermonde matrix

// Test 1.1.3: Mass matrix properties
func TestMassMatrixProperties(t *testing.T)
// - Compute M_ij = <ψ_i, ψ_j>
// - Verify M = I (identity) for orthonormal basis
// - Check sparsity pattern for nodal basis
```

### Step 1.2: Derivative Matrices
**Implementation:**
```go
// Derivative operators in collapsed coordinates
func (tb *TetBasis) ComputeGradientMatrices() (Dr, Ds, Dt [][]float64)

// Transform derivatives to physical coordinates
func (tb *TetBasis) PhysicalDerivatives(Dr, Ds, Dt [][]float64, J [3][3]float64) (Dx, Dy, Dz [][]float64)

// Strong form differentiation matrices
func (tb *TetBasis) StrongFormOperators() (Dr, Ds, Dt [][]float64)
```

**Unit Tests:**
```go
// Test 1.2.1: Polynomial derivative exactness
func TestTetGradientAccuracy(t *testing.T)
// - Test ∂/∂x (x^i*y^j*z^k) for i+j+k <= P-1
// - Verify exactness at interpolation nodes
// - Maximum error < 1e-14

// Test 1.2.2: Divergence operator accuracy
func TestTetDivergenceAccuracy(t *testing.T)
// - Test ∇·F for polynomial vector fields F
// - Verify discrete divergence theorem
// - Check conservation properties
```

### Step 1.3: Mode Ordering and Access
**Implementation:**
```go
// Mode ordering by polynomial degree (hierarchical)
type ModeOrdering struct {
    ModeToIJK map[int][3]int // mode index -> (i,j,k)
    IJKToMode map[[3]int]int // (i,j,k) -> mode index
}

func (tb *TetBasis) GenerateModeOrdering() *ModeOrdering
func (tb *TetBasis) GetModesOfOrder(p int) []int // Return indices of modes with total order p
func (tb *TetBasis) GetModalSlice(pMin, pMax int) []int // Modes with order in [pMin, pMax]
```

**Unit Tests:**
```go
// Test 1.3.1: Mode ordering consistency
func TestTetModeOrdering(t *testing.T)
// - Verify bijection between mode index and (i,j,k)
// - Check hierarchical ordering: lower orders come first
// - Verify mode count per order: (p+1)(p+2)/2 for order p

// Test 1.3.2: Modal truncation
func TestModalTruncation(t *testing.T)
// - Project high-order function to lower order
// - Verify optimal L2 approximation property
// - Test spectral convergence rates
```

### Step 1.4: Warburton Interpolation Nodes
**Implementation:**
```go
// Generate optimized interpolation nodes
// Reference: Warburton (2006), "An explicit construction of interpolation nodes on the simplex"
func (tb *TetBasis) GenerateWarpBlendNodes() (r, s, t []float64)
func WarpFunction(n int, r []float64) []float64
func BlendFunction(a, b, c float64) float64
```

**Unit Tests:**
```go
// Test 1.4.1: Node distribution quality
func TestWarpBlendNodes(t *testing.T)
// - Compute Lebesgue constant
// - Verify edge nodes match 1D Gauss-Lobatto
// - Check face nodes match 2D warp-blend triangle nodes

// Test 1.4.2: Interpolation stability
func TestInterpolationStability(t *testing.T)
// - Test Runge function interpolation
// - Verify spectral convergence
// - Compare condition numbers with equispaced nodes
```

---

## Phase 2: Hexahedral Elements (Tensor Product Basis with Line-Based Approach)

### Step 2.1: Tensor Product Lagrange Basis
**Implementation:**
```go
type HexBasis struct {
    Order     int
    Np        int // = (P+1)^3
    Nodes1D   []float64 // Gauss-Lobatto-Legendre nodes
    Weights1D []float64 // GLL quadrature weights
}

// Generate GLL nodes and weights
func GenerateGLLNodesWeights(order int) (nodes, weights []float64)

// Lagrange basis on GLL nodes
func LagrangeBasis1D(i int, x float64, nodes []float64) float64

// Tensor product basis: ψ_ijk(r,s,t) = L_i(r) * L_j(s) * L_k(t)
func (hb *HexBasis) EvaluateTensorBasis(r, s, t float64) []float64
```

**Unit Tests:**
```go
// Test 2.1.1: GLL properties
func TestGLLNodesWeights(t *testing.T)
// - Verify nodes include endpoints [-1,1]
// - Check quadrature exactness for polynomials up to 2P-1
// - Verify symmetry of nodes and weights

// Test 2.1.2: Tensor basis properties  
func TestHexTensorBasis(t *testing.T)
// - Verify Kronecker delta property at nodes
// - Check partition of unity
// - Test tensor product structure
```

### Step 2.2: Line-Based DG Implementation
**Implementation:**
```go
// Line-based operators following Persson (2013)
// Reference: Persson (2013), "A Sparse and High-Order Accurate Line-Based 
// Discontinuous Galerkin Method for Unstructured Meshes", J. Comput. Phys. 233, 414-429
type LineDGOperator struct {
    Order int
    D1D   [][]float64 // 1D derivative matrix
    M1D   [][]float64 // 1D mass matrix
    F1D   [][]float64 // 1D flux matrix
}

// Apply derivative along single direction (sparse operation)
func (hb *HexBasis) ApplyLineDerivative(u []float64, direction int) []float64

// Line-based flux integral with O(P) connectivity per node
func (hb *HexBasis) LineBasedFluxIntegral(flux []float64, normal [3]float64) []float64

// Generate sparse connectivity pattern
func (hb *HexBasis) ComputeLineConnectivity() [][]int
```

**Unit Tests:**
```go
// Test 2.2.1: Sparsity verification
func TestLineBasedSparsity(t *testing.T)
// - Count non-zeros per row: should be O(P) not O(P²)
// - Verify line-only connectivity pattern
// - Compare memory usage with standard DG

// Test 2.2.2: Accuracy comparison
func TestLineBasedAccuracy(t *testing.T)
// - Compare with standard DG on test problems
// - Verify identical results (up to round-off)
// - Test conservation properties
```

### Step 2.3: Sum-Factorization Implementation
**Implementation:**
```go
// Efficient tensor operations via sum-factorization
// Complexity: O(P^4) in 3D instead of O(P^6)
func (hb *HexBasis) SumFactorizedGradient(u []float64) (ux, uy, uz []float64)

// Tensor contraction: apply 1D operator in one direction
func TensorContraction3D(u []float64, op1D [][]float64, direction int, P int) []float64

// Kronecker product representation
func (hb *HexBasis) DecomposeOperator(M [][]float64) (Mr, Ms, Mt [][]float64)
```

**Unit Tests:**
```go
// Test 2.3.1: Sum-factorization correctness
func TestSumFactorization(t *testing.T)
// - Compare with direct matrix-vector product
// - Verify to machine precision
// - Measure FLOP count reduction

// Test 2.3.2: Performance scaling
func BenchmarkSumFactorization(b *testing.B)
// - Benchmark vs direct multiplication
// - Verify O(P^4) vs O(P^6) scaling
// - Test different polynomial orders
```

### Step 2.4: Tensor Product Preconditioning
**Implementation:**
```go
// SVD-based tensor approximation following Pazner & Persson (2018)
// Reference: Pazner & Persson (2018), "Approximate tensor-product preconditioners 
// for very high order discontinuous Galerkin methods", J. Comput. Phys. 354, 344-369
type TensorPreconditioner struct {
    Order      int
    Tolerance  float64
    Factors    []KroneckerFactor
}

type KroneckerFactor struct {
    U1, U2, U3 [][]float64
    Sigma      []float64
}

// Approximate block Jacobi via SVD
func ApproximateTensorProduct(J [][]float64, tol float64) *TensorPreconditioner

// Apply preconditioner: O(P^2) storage, O(P^4) work
func (tp *TensorPreconditioner) Apply(r []float64) []float64
```

**Unit Tests:**
```go
// Test 2.4.1: Approximation quality
func TestTensorApproximation(t *testing.T)
// - Verify ||J - J_approx||_F / ||J||_F < tolerance
// - Check storage reduction from O(P^6) to O(P^2)
// - Test various tolerance levels

// Test 2.4.2: Preconditioning effectiveness
func TestTensorPreconditioning(t *testing.T)
// - Compare GMRES iterations with exact block Jacobi
// - Test on Poisson and advection-diffusion
// - Verify scalability with P
```

### Step 2.5: Mode Ordering for Hex Elements
**Implementation:**
```go
// Lexicographic mode ordering for tensor products
func (hb *HexBasis) GenerateTensorModeOrdering() *ModeOrdering

// Access modes by directional orders
func (hb *HexBasis) GetModesWithOrders(pr, ps, pt int) []int

// Modal filtering by direction
func (hb *HexBasis) ApplyDirectionalFilter(u []float64, filterR, filterS, filterT func(int) float64) []float64
```

**Unit Tests:**
```go
// Test 2.5.1: Tensor mode structure
func TestHexModeOrdering(t *testing.T)
// - Verify tensor product ordering
// - Check mode count: (pr+1)(ps+1)(pt+1)
// - Test directional mode access

// Test 2.5.2: Directional filtering
func TestDirectionalFiltering(t *testing.T)
// - Apply different filters per direction
// - Verify selective attenuation
// - Test on discontinuous data
```

---

## Phase 3: Prismatic Elements (Hybrid Tensor Product)

### Step 3.1: Triangle-Line Tensor Product Basis
**Implementation:**
```go
type PrismBasis struct {
    TriangleOrder int
    LineOrder     int  
    Np            int // = (P+1)(P+2)/2 * (Q+1)
    TriBasis      *TriangleBasis // 2D PKD basis
    LineBasis     *LineBasis1D    // 1D GLL Lagrange basis
}

// Basis functions: ψ_ijk(r,s,t) = φ^tri_ij(r,s) * L_k(t)
func (pb *PrismBasis) EvaluateBasis(r, s, t float64) []float64

// Collapsed coordinates for triangle
func TriangleCollapse(r, s float64) (a, b float64) {
    // a = 2(1+r)/(1-s) - 1 if s != 1, else a = -1
    // b = s
}
```

**Unit Tests:**
```go
// Test 3.1.1: Hybrid basis structure
func TestPrismHybridBasis(t *testing.T)
// - Verify tensor product form
// - Check total DOF count
// - Test basis completeness

// Test 3.1.2: Mixed orthogonality
func TestPrismOrthogonality(t *testing.T)
// - Verify orthogonality in z-direction (GLL)
// - Check PKD orthogonality in triangle
// - Test mass matrix block structure
```

### Step 3.2: Efficient Line Operations
**Implementation:**
```go
// Line operations in vertical direction
func (pb *PrismBasis) VerticalDerivative(u []float64) []float64

// Exploit tensor structure for z-integration
func (pb *PrismBasis) EfficientZQuadrature(f []float64) []float64

// Separate vertical modes for efficient operations
func (pb *PrismBasis) SeparateVerticalModes(u []float64) [][]float64
```

**Unit Tests:**
```go
// Test 3.2.1: Vertical efficiency
func TestPrismVerticalOps(t *testing.T)
// - Verify O(PQ) complexity for z-derivatives
// - Compare timing with full 3D operations
// - Check accuracy preservation

// Test 3.2.2: Mode separation
func TestPrismModeSeparation(t *testing.T)
// - Test separation/recombination identity
// - Verify no information loss
// - Check efficiency gains
```

### Step 3.3: Mode Ordering for Prisms
**Implementation:**
```go
// Hierarchical ordering: triangle modes tensored with line modes
type PrismModeOrdering struct {
    TriangleModes []int
    LineModes     []int
    PrismToTriLine map[int][2]int // prism mode -> (tri_mode, line_mode)
}

func (pb *PrismBasis) GeneratePrismModeOrdering() *PrismModeOrdering

// Get all modes with specific triangle order p and line order q
func (pb *PrismBasis) GetModesOfOrders(pTri, qLine int) []int
```

**Unit Tests:**
```go
// Test 3.3.1: Prism mode structure
func TestPrismModeOrdering(t *testing.T)
// - Verify hierarchical ordering
// - Check mode separation by orders
// - Test reconstruction from components

// Test 3.3.2: Anisotropic truncation
func TestAnisotropicTruncation(t *testing.T)
// - Truncate different orders in triangle/line
// - Verify optimal approximation
// - Test directional convergence rates
```

### Step 3.4: Tensor Preconditioning for Prisms
**Implementation:**
```go
// Exploit separable structure when possible
type PrismPreconditioner struct {
    TrianglePrecond *TrianglePreconditioner  
    LinePrecond     *LinePreconditioner1D
    CouplingMatrix  *SparseCoupling // For non-separable terms
}

func BuildPrismPreconditioner(J [][]float64) *PrismPreconditioner
```

**Unit Tests:**
```go
// Test 3.4.1: Separable problems
func TestPrismSeparablePrecond(t *testing.T)
// - Test on Laplacian (separable)
// - Verify convergence improvement
// - Check memory efficiency

// Test 3.4.2: Coupled problems  
func TestPrismCoupledPrecond(t *testing.T)
// - Test on advection-diffusion
// - Compare with full preconditioner
// - Verify robustness
```

---

## Phase 4: Pyramidal Elements

### Step 4.1: Rational Basis Functions
**Implementation:**
```go
type PyramidBasis struct {
    Order int
    Np    int // = (P+1)(P+2)(2P+3)/6
}

// Modified coordinates to handle apex singularity
func PyramidCollapse(r, s, t float64) (xi, eta, zeta float64) {
    // xi = r*2/(1-t) if t != 1, else 0
    // eta = s*2/(1-t) if t != 1, else 0
    // zeta = t
}

// Rational basis with proper normalization
func (pyb *PyramidBasis) EvaluateRationalBasis(r, s, t float64) []float64
```

**Unit Tests:**
```go
// Test 4.1.1: Apex handling
func TestPyramidApexBehavior(t *testing.T)
// - Evaluate basis near apex (t → 1)
// - Verify finite values
// - Check continuity

// Test 4.1.2: Basis properties
func TestPyramidBasisProperties(t *testing.T)
// - Verify polynomial space dimension
// - Test orthogonality where applicable
// - Check conditioning away from apex
```

### Step 4.2: Special Quadrature
**Implementation:**
```go
// Quadrature avoiding apex singularity
func (pyb *PyramidBasis) GenerateQuadratureRule() (r, s, t, w []float64)

// Integration with singularity handling
func (pyb *PyramidBasis) IntegrateFunction(f func(r, s, t float64) float64) float64
```

**Unit Tests:**
```go
// Test 4.2.1: Quadrature accuracy
func TestPyramidQuadrature(t *testing.T)
// - Integrate polynomials exactly
// - Test near-apex behavior
// - Verify positive weights

// Test 4.2.2: Mass matrix conditioning
func TestPyramidMassMatrix(t *testing.T)
// - Compute condition number
// - Verify positive definiteness
// - Compare with tet/hex elements
```

### Step 4.3: Mode Ordering for Pyramids
**Implementation:**
```go
// Special ordering due to rational basis
func (pyb *PyramidBasis) GeneratePyramidModeOrdering() *ModeOrdering

// Access modes avoiding apex issues
func (pyb *PyramidBasis) GetStableModes(maxOrder int) []int
```

**Unit Tests:**
```go
// Test 4.3.1: Pyramid mode structure
func TestPyramidModeOrdering(t *testing.T)
// - Verify mode count formula
// - Check hierarchical property
// - Test stability of high modes

// Test 4.3.2: Modal filtering
func TestPyramidModalFilter(t *testing.T)
// - Apply filtering near apex
// - Verify stability improvement
// - Test accuracy preservation
```

---

## Key References

1. **Hesthaven, J.S. and Warburton, T.** (2008). "Nodal Discontinuous Galerkin Methods: Algorithms, Analysis, and Applications", Springer. *[Primary reference for PKD basis and nodal DG methods]*

2. **Persson, P.-O.** (2013). "A Sparse and High-Order Accurate Line-Based Discontinuous Galerkin Method for Unstructured Meshes", J. Comput. Phys., 233, 414-429. *[Line-based DG approach reducing connectivity]*

3. **Pazner, W. and Persson, P.-O.** (2018). "Approximate tensor-product preconditioners for very high order discontinuous Galerkin methods", J. Comput. Phys., 354, 344-369. *[SVD-based tensor preconditioning]*

4. **Warburton, T.** (2006). "An explicit construction of interpolation nodes on the simplex", J. Eng. Math., 56, 247-262. *[Warp-blend nodes for simplices]*

5. **Fortunato, D., Hale, N., and Townsend, A.** (2021). "The ultraspherical spectral element method", J. Comput. Phys., 436, 110087. *[Modern spectral element approaches]*
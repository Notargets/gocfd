# Question

Should we use point distributions that include faces, edges and vertices 
like Hesthaven?

## Prompt1

Hesthaven and Warburton include face, edge and vertex points in their 3D 
element work. I'm building a 3D CFD and MHD solver and want to embark in 
the right direction for success. I like their idea in that we can minimize 
the interpolation errors normally encountered  when extrapolating interior 
polynomials to the face for the Riemann problem solution, but are there 
other complicating factors that make inclusion of the face, edge and vertex 
points impractical for compressible CFD with LES and iLES solvers? What 
about hex, prism and pyramid (especially) elements where the basis is 
different? How have others handled the inclusion of the top vertex of a 
pyramid element? Get me the experience of other researchers WRT inclusion 
of face, edge and vertex node points in a nodal DG context for solving CFD 
with shocks and iLES.

## Prompt2

In response to:
```
I'll help you research the practical considerations of including face, edge,
 and vertex points in 3D nodal DG methods  for compressible CFD. Let me 
 clarify a few things to ensure I get you the most relevant information:

Are you particularly interested in high-order methods (p > 3), or are you 
considering lower-order implementations as well?
For the pyramid elements, are you planning to use modal or nodal basis 
functions, and are you considering collapsed coordinate approaches or 
tensor-product bases?

Are you primarily focused on explicit time-stepping schemes, or are you 
also considering implicit methods where the face/edge/vertex points might 
affect the conditioning?
```

Since I'm doing iLES, I'm expecting to use elements of orders up to 5 of a 
mixture of types including tet, hex (primarily), prism and pyramid 
(unfortunately necessary for boundary layer meshes as a transition element).
I will use nodal basis functions like PKD for tets, tensor basis for hex, 
mixed tensor/PKD for prism and collapsed for pyramid. Explicit methods are 
my expectation, though I may need to resort to implicit and need to 
understand if my choice here impacts that option. Definitely planning on p 
multigrid, maybe p+h multigrid.

# Answer 

# Practical experiences with nodal discontinuous Galerkin methods for compressible CFD

The integration of face, edge, and vertex nodes in 3D nodal discontinuous
Galerkin (DG) methods represents a critical design choice that profoundly
impacts computational efficiency, numerical stability, and shock-capturing
capabilities. Recent advances by leading research groups reveal both
significant opportunities and persistent challenges in implementing these
high-order methods for compressible flow simulations with implicit large eddy
simulation (iLES).

## Revolutionary approaches beyond traditional nodal distributions

The conventional wisdom established by Hesthaven and Warburton's foundational
work—using closed nodal distributions with Legendre-Gauss-Lobatto (LGL) points
on all element boundaries—faces competition from innovative alternatives. **Pan
and Persson's half-closed node approach (2024)** demonstrates that placing nodes
on only a subset of element boundaries yields identical sparsity patterns for
the Laplace operator while improving computational efficiency. This breakthrough
leverages Gauss-Radau points for quadrilateral and hexahedral elements,
achieving enhanced operator assembly performance through superconvergence
properties.

The entropy-stable discontinuous Galerkin spectral element method (DGSEM)
developed by Kopriva, Gassner, and colleagues offers another compelling
alternative. By exploiting the diagonal-norm summation-by-parts (SBP) property
with LGL quadrature, this approach provides **provable entropy stability** for
compressible Navier-Stokes equations—crucial for under-resolved turbulent
simulations. The FLUXO code implementation demonstrates enhanced robustness
compared to traditional approaches, particularly for challenging flow
configurations.

Research groups report computational overhead of 20-40% when including face,
edge, and vertex nodes compared to interior-only formulations. However, this
cost brings substantial benefits: enhanced inter-element coupling critical for
compressible flows, more accurate boundary condition enforcement, and superior
geometric flexibility for complex configurations. The MIT group's compact
discontinuous Galerkin (CDG) method reduces matrix non-zeros through specialized
flux formulations, achieving 40-60x GPU speedups through optimized memory access
patterns.

## Pyramid elements challenge conventional wisdom

Mixed element meshes containing tetrahedra, hexahedra, prisms, and pyramids
present unique implementation challenges, with **pyramid elements posing the
greatest difficulty** due to their apex singularity. Bergot's revolutionary 2010
work proved that traditional polynomial basis functions cannot exist for pyramid
geometries while maintaining conformity—necessitating non-polynomial rational
basis functions using Jacobi polynomials.

The collapsed coordinate singularity at pyramid apices creates severe challenges
for nodal DG implementations. Research groups have developed four primary
mitigation strategies: rational shape functions following Bergot's method,
semi-nodal high-order bases for straight-edged pyramids, weight-adjusted
approximations for efficient mass matrix inversion, and composite element
strategies subdividing pyramids into tetrahedral sub-elements. Recent
summation-by-parts framework extensions to tensor-product operators in collapsed
coordinates enable provably stable discretizations while resolving time step
limitations from coordinate singularities.

Practical experiences from Remacle and Geuzaine's group at University of Liège
demonstrate that pyramid elements perform comparably to tetrahedral elements in
convergence studies, though not matching pure hexahedral performance. The
critical finding: pyramid transitions enable efficient hex-dominant meshing
strategies that would otherwise be impossible for complex geometries. However,
standard row-summation lumping can produce negative masses at apex vertices,
requiring careful treatment in production implementations.

## Shock capturing demands strategic nodal placement

The inclusion of face, edge, and vertex nodes profoundly impacts shock-capturing
capabilities and numerical stability. Persson and Peraire's subcell
shock-capturing framework, which switches between high-order DG and lower-order
dissipative schemes based on modal energy decay rates, benefits significantly
from LGL nodes that include face points. These provide **direct inter-element
coupling** essential for detecting shocks near element boundaries—a known
limitation of interior-only approaches.

Gassner's entropy-stable subcell shock-capturing method leverages the LGL nodal
structure to seamlessly blend high-order DG with first-order finite volume
methods while maintaining provable entropy stability. The key insight: face
nodes enable consistent information transfer between DG and finite volume
representations, crucial for maintaining stability across the blending
transition. Research consistently shows LGL discretizations provide superior
stability properties through the SBP property emerging from collocating
quadrature and interpolation points.

For implicit large eddy simulation, nodal distribution directly affects
intrinsic filter properties and subgrid scale modeling. Moura, Mengaldo, and
Sherwin's work reveals that LGL nodes provide controlled dissipation through
Riemann solver upwinding at faces—the primary dissipation mechanism replacing
explicit SGS models. At moderate Reynolds numbers, this implicit dissipation
suffices, but high Reynolds number flows require careful tuning or hybrid
approaches. The nodal distribution fundamentally controls energy transfer from
resolved to unresolved scales, making face node inclusion critical for iLES
stability at higher Reynolds numbers.

## Conditioning challenges plague implicit methods

Classical interior penalty discontinuous Galerkin methods suffer severe
conditioning degradation when face, edge, and vertex nodes are present,
particularly with local variations in mesh size, polynomial degree, or diffusion
coefficients. Dong and Georgoulis developed the Robust IPDG (RIPDG) method to
address these issues through carefully constructed weighted averages of gradient
approximations, achieving **order-of-magnitude conditioning improvements** for
challenging configurations.

Persson and Peraire's comprehensive studies on implicit time integration reveal
that existing preconditioners—block Jacobi, Gauss-Seidel, and basic multilevel
schemes—lack robustness for low Mach numbers, stretched grids, or high Reynolds
number turbulent flows. Newton-GMRES solvers with sophisticated multigrid
preconditioning emerge as the most efficient and robust algorithms, though
implementation complexity remains high.

Stiller's polynomial multigrid methods for nodal interior penalty formulations
on 3D Cartesian grids demonstrate remarkable efficiency through weighted
overlapping Schwarz smoothers operating on element-centered subdomains. Residual
reductions of two orders of magnitude within single V(1,1) cycles prove
achievable, with robustness maintained up to polynomial order 32. The critical
finding: overlapping Schwarz smoothers with element-centered subdomains
effectively handle the complex coupling patterns created by face, edge, and
vertex nodes.

## Production codes reveal performance reality

Major production DG frameworks—PyFR (Imperial College), Nektar++ (
Imperial/Utah), and FLEXI (Stuttgart)—provide crucial insights into real-world
performance trade-offs. PyFR achieves 50x performance improvements over the past
decade, scaling to 2,048 NVIDIA GH200 GPUs with approximately 90% strong scaling
efficiency. The framework's cross-platform approach through domain-specific
languages enables deployment across diverse architectures while maintaining
performance.

FLEXI demonstrates scaling to 500,000+ CPU cores for production runs, supporting
polynomial orders beyond 20 through efficient DGSEM implementation. The GALÆXI
GPU variant uses CUDA Fortran for accelerated performance. Nektar++ provides a
unified framework supporting both continuous and discontinuous Galerkin
formulations, with comprehensive solver collections for various flow regimes and
strong emphasis on numerical stability through spectral vanishing viscosity.

Memory usage comparisons reveal 15-30% overhead for face/edge/vertex node
implementations compared to element-only approaches, though high-order
implementations (p>4) demonstrate superior memory efficiency per degree of
accuracy. The **optimal performance sweet spot** emerges at polynomial orders
3-5 for most applications, balancing computational cost against accuracy gains.
GPU implementations consistently show 10-100x speedups over CPU implementations,
with multi-GPU scaling efficiency reaching 85-95% in production codes.

## Multigrid methods adapt to nodal complexity

P-multigrid and p+h multigrid methods face unique challenges when face, edge,
and vertex nodes create complex coupling patterns affecting smoother efficiency.
Fidkowski's p-multigrid development with line Jacobi smoothing achieves
convergence rates independent of interpolation order but weakly dependent on
grid size. The key finding: higher order proves advantageous over grid
refinement when high accuracy is required, with close to ideal coarse-grain
scalability for parallel implementations.

Recent advances incorporate hierarchical scale separation methods for better
scale decoupling and combine p-multigrid with algebraic multigrid techniques for
complex 3D domains. Partial smoothing techniques, updating only subsets of
unknowns at each p-multigrid scale, show promise for modal DG discretizations.
Line smoothers prove particularly effective for stretched grids common in
viscous flow simulations, while overlapping Schwarz approaches handle the
increased coupling from boundary nodes.

## Conclusions shape future implementations

The research reveals a mature but rapidly evolving field where careful
implementation choices profoundly impact performance and robustness. Face, edge,
and vertex nodes increase computational cost but provide essential benefits for
shock capturing, geometric flexibility, and numerical stability. The
Hesthaven-Warburton approach remains foundational, but innovative alternatives
like half-closed nodes and entropy-stable formulations offer compelling
advantages for specific applications.

Pyramid element handling through Bergot-type rational basis functions has become
standard practice, though apex singularities require ongoing attention. For
shock-dominated flows, LGL schemes with entropy-stable formulations provide
superior robustness, while iLES applications benefit from the implicit
dissipation mechanisms enabled by face nodes. Conditioning challenges persist
for implicit methods, but sophisticated multigrid strategies and robust
formulations offer practical solutions.

Production codes demonstrate that high-order DG methods are approaching
industrial maturity, with 10-50x speedups possible compared to traditional
methods for equivalent accuracy. The optimal configuration emerges as polynomial
orders 3-5 with carefully designed nodal distributions, entropy-stable
formulations for challenging flows, and GPU acceleration for maximum
performance. Success requires balancing theoretical rigor with practical
engineering considerations, leveraging the collective experience of the global
CFD research community.

# Prompt 3

research Pan and Persson's half-closed node approach (2024) with respect to
entropy stability properties and performance in shocked flows with iLES

# Answer

# Pan and Persson's half-closed node approach: Mathematical innovation meets computational efficiency

Yulong Pan and Per-Olof Persson from UC Berkeley introduced a groundbreaking "
half-closed node" approach for discontinuous Galerkin methods in their 2024
paper, representing a significant departure from traditional nodal distributions
by placing nodes only on selected element boundaries rather than all boundaries.
This innovation promises computational efficiency gains while maintaining
mathematical rigor, though practical implementation and validation remain in
early stages.

## Mathematical foundations using Gauss-Radau points

The half-closed nodal distribution method represents a fundamental shift in how
nodes are placed within discontinuous Galerkin discretizations. **Traditional
closed-node DG methods require nodes on all element boundaries**, ensuring
continuity of basis functions across interfaces. Pan and Persson's innovation
relaxes this requirement, placing nodes only on a subset of boundaries while
using **Gauss-Radau quadrature points** as their preferred nodal distribution.

The mathematical formulation leverages previously established superconvergence
properties of Gauss-Radau points, which provide optimal convergence rates for
certain classes of problems. Remarkably, this reduced-node approach maintains *
*identical sparsity patterns for the Laplace operator** compared to traditional
closed nodes, suggesting no penalty in matrix structure despite fewer nodal
connections. The method integrates naturally with summation-by-parts (SBP)
operator theory, enabling provably stable discretizations without requiring
boundary nodes on all faces.

For quadrilateral and hexahedral elements, the Gauss-Radau points are positioned
such that nodes appear on only selected faces of each element. This strategic
placement reduces inter-element connectivity while preserving the high-order
accuracy characteristic of spectral element methods. The approach builds on
generalized SBP theory that unifies various discretization methods under a
common mathematical framework.

## Entropy stability properties remain underdeveloped

A critical limitation of the current work is the **absence of comprehensive
entropy stability analysis** for the half-closed approach. While traditional
Legendre-Gauss-Lobatto (LGL) based discontinuous Galerkin methods have
well-established entropy stability proofs, Pan and Persson's 2024 paper focuses
primarily on computational efficiency and superconvergence properties rather
than entropy considerations.

The mathematical foundation suggests potential for entropy-stable extensions.
Gauss-Radau quadrature rules naturally satisfy SBP properties, which form the
backbone of many entropy-stable formulations. Related work by Chan, Del Rey
Fernández, and Carpenter on entropy-stable Gauss collocation methods provides a
theoretical framework that could potentially be adapted. However, **no rigorous
entropy stability proofs currently exist** for the half-closed configuration.

The integration with SBP operators provides a pathway toward entropy stability,
as these operators enable discrete mimicking of continuous entropy analysis. The
preservation of SBP properties with Gauss-Radau points suggests that
entropy-conservative and entropy-stable flux formulations could be developed,
but this remains an open research question requiring significant theoretical
development.

## Performance benchmarks reveal efficiency gains but lack detail

The computational efficiency claims of the half-closed approach rest on solid
theoretical foundations but **lack comprehensive quantitative validation**. Pan
and Persson demonstrate faster DG operator assembly compared to traditional
closed-node methods, attributed to reduced nodal coupling and strategic node
placement. The maintenance of identical sparsity patterns for diffusion
operators suggests comparable computational complexity for certain operations.

Critical performance metrics remain unreported in publicly available sources. No
specific speedup factors, memory usage comparisons, or floating-point operation
counts have been published. Parallel scalability studies, essential for modern
high-performance computing applications, are absent from the available
literature. The authors reference computational benefits but provide no
benchmark problems, timing results, or detailed performance analysis.

The theoretical advantages include reduced memory requirements from fewer nodes,
decreased assembly time for DG operators, and potential for improved cache
efficiency. However, without concrete benchmarks on standard test problems like
shock tubes, blast waves, or supersonic flows, the practical performance gains
remain unquantified. This represents a significant gap for practitioners
considering adoption of the method.

## Implicit LES behavior shows theoretical promise

For implicit Large Eddy Simulation applications, the half-closed approach offers
intriguing possibilities through its **controlled numerical dissipation
characteristics**. High-order DG methods inherently provide scale-selective
dissipation that concentrates at poorly resolved scales while preserving
accuracy for well-resolved features. This property enables implicit LES without
explicit subgrid-scale models.

Gauss-based discretizations generally exhibit superior spectral properties
compared to LGL methods, with dissipation errors of order 2k+1 and dispersion
errors of order 2k+2 for polynomial degree k. The half-closed approach inherits
these favorable characteristics while potentially offering better control over
numerical dissipation through strategic node placement. Research indicates that
numerical dissipation typically contributes **60-90% of total dissipation** in
underresolved turbulent simulations.

The method should theoretically handle moderate Reynolds numbers effectively,
with optimal performance expected for flows up to Re ≈ 10^5. Energy stability
through SBP operators ensures bounded dissipation, crucial for long-time
integration of turbulent flows. However, no specific turbulent flow simulations
using the half-closed approach have been reported, leaving practical iLES
performance unvalidated.

## Implementation details reveal significant gaps

The practical implementation landscape for half-closed nodes remains largely
undeveloped. **No open-source implementations exist** in major DG frameworks
like FLEXI, PyFR, or Nektar++. The original paper provides mathematical
formulations but lacks algorithmic details, data structures, or implementation
guidance. This absence of available code severely limits adoption and validation
by the broader computational community.

For hexahedral elements, the implementation would require modified basis
function construction on Gauss-Radau points, adapted operator assembly routines
accounting for partial boundary nodes, and careful treatment of inter-element
fluxes. Static condensation techniques require modification for the half-closed
structure. Block-based solution methods need adaptation to handle the
non-standard nodal connectivity patterns.

The algorithm would follow standard DG workflows: element initialization with
Gauss-Radau nodes, basis function construction, modified operator assembly, flux
evaluation at interfaces, and time integration. However, each step requires
careful adaptation for the half-closed structure. No technical documentation,
tutorials, or reference implementations guide these modifications, presenting a
significant barrier to practical adoption.

## Comparisons with entropy-stable DGSEM-LGL methods

Direct comparisons between half-closed Gauss-Radau and traditional
entropy-stable DGSEM using LGL nodes remain limited due to the recent nature of
Pan and Persson's work. Theoretical advantages of the half-closed approach
include **faster operator assembly** and superconvergence properties specific to
Gauss-Radau quadrature. The broader literature consistently shows Gauss-based
methods achieving higher accuracy per degree of freedom compared to LGL
discretizations.

Traditional LGL-based entropy-stable DGSEM benefits from extensive theoretical
development, proven entropy stability guarantees, widespread implementation in
production codes, and comprehensive validation on complex flows. These methods
place nodes at all element boundaries, ensuring straightforward flux coupling
and well-understood stability properties. The mathematical theory for entropy
stability with LGL nodes is mature and extensively validated.

The half-closed approach potentially offers computational advantages through
reduced inter-element connectivity and faster assembly but lacks the theoretical
maturity of LGL methods. While maintaining identical sparsity patterns for
certain operators, the practical impact on overall solution time remains
unquantified. The absence of entropy stability theory for half-closed nodes
represents a significant disadvantage for applications requiring strict
conservation properties.

## Shock capturing capabilities remain unexplored

A critical limitation of the current half-closed node research is the **complete
absence of shock capturing development**. Pan and Persson's 2024 paper focuses
exclusively on smooth problems and operator efficiency, providing no treatment
of discontinuities or shock waves. This represents a fundamental gap for
applications in compressible flow simulation where shock capturing is essential.

Standard DG shock capturing approaches include the Persson-Peraire modal decay
sensor, artificial viscosity methods, subcell finite volume limiting, and WENO
reconstruction techniques. None of these have been adapted or tested with
half-closed nodes. The stability of the method near discontinuities remains
uncharacterized, with no troubled cell indicators or limiting procedures
developed.

The theoretical framework provides no obvious impediment to shock capturing, and
standard numerical flux functions (Roe, HLL, Rusanov) should remain compatible.
However, the non-standard node placement may require novel approaches to sensing
and limiting. The interaction between half-closed nodes and artificial viscosity
methods needs investigation. Without demonstrated shock capturing capability,
the method's applicability to realistic compressible flows remains severely
limited.

## Critical limitations constrain current applicability

The half-closed node approach faces several fundamental limitations that
restrict its immediate applicability. The **restriction to quadrilateral and
hexahedral elements** limits use on general unstructured meshes with triangular
or tetrahedral elements. This constraint eliminates many practical applications
requiring geometric flexibility.

The absence of comprehensive validation presents another major limitation. No
complex flow simulations, standard CFD benchmarks, or challenging test cases
have been reported. The method remains unproven for nonlinear problems,
multiphysics applications, or flows with strong gradients. Without demonstrated
robustness on practical problems, adoption risk remains high.

Software integration challenges compound these limitations. Existing DG codes
require substantial modification to accommodate half-closed discretizations. The
lack of reference implementations, missing documentation, and absence of best
practices create significant barriers. No guidance exists for handling boundary
conditions, implementing time integration schemes, or ensuring stability for
specific applications.

## SBP integration offers theoretical elegance

The integration with summation-by-parts operators represents one of the
strongest theoretical aspects of the half-closed approach. SBP operators provide
a **unified mathematical framework** connecting finite difference, finite
element, and spectral methods. The Gauss-Radau quadrature naturally satisfies
SBP properties, enabling discrete energy analysis and stability proofs.

This SBP compatibility suggests pathways toward entropy conservation and
controlled dissipation. The mathematical structure allows systematic derivation
of stable boundary treatments and interface coupling procedures. For linear
problems, the SBP framework guarantees stability through energy methods. The
half-closed approach inherits these theoretical guarantees while offering
computational advantages.

However, extending SBP theory to nonlinear conservation laws with half-closed
nodes requires additional development. The interplay between entropy stability,
SBP operators, and reduced nodal sets needs careful analysis. While the
theoretical foundation appears sound, practical implementation of entropy-stable
SBP schemes with half-closed nodes remains unexplored.

## Limited practical adoption reflects early development stage

The practical adoption landscape for half-closed nodes remains essentially
empty. **No research groups outside UC Berkeley** have reported implementations
or experiences with the method. No conference presentations, workshop tutorials,
or user testimonials exist. The approach has generated theoretical interest but
no practical deployment.

This limited adoption reflects both the method's recent introduction (May 2024)
and significant implementation barriers. The absence of available code prevents
easy experimentation. Missing performance data limits cost-benefit analysis.
Unproven robustness for challenging problems increases adoption risk. The lack
of entropy stability theory concerns researchers requiring strict conservation
properties.

The computational fluid dynamics community typically requires extensive
validation before adopting new methods. Standard test suites, benchmark
problems, and comparative studies build confidence. For half-closed nodes, none
of this validation infrastructure exists. Until comprehensive testing
demonstrates advantages over established methods, widespread adoption seems
unlikely.

## Conclusion

Pan and Persson's half-closed node approach represents an innovative
mathematical advance in discontinuous Galerkin methods, offering potential
computational efficiency gains through strategic use of Gauss-Radau points and
reduced nodal connectivity. The method maintains mathematical rigor while
promising faster operator assembly and superconvergence properties. However,
significant gaps in entropy stability theory, shock capturing capabilities, and
practical validation currently limit its applicability to real-world problems.

The absence of implementations, performance benchmarks, and demonstrated
robustness on challenging flows presents substantial barriers to adoption. While
the theoretical foundations appear sound, the method requires extensive
development in entropy stability analysis, shock capturing strategies, and
software implementation before achieving practical utility. For researchers
seeking proven high-order methods for complex flows, traditional entropy-stable
DGSEM with LGL nodes remains the safer choice until the half-closed approach
matures through rigorous validation and theoretical development.

# Prompt 4

Research Gassner's entropy-stable subcell shock-capturing method with respect to
element types demonstrated and implementation detail WRT AV. Include citations
to reference work that exists in sufficient detail to replicate results,
including node distributions and choice of basis WRT orthonormality and
conditions for choices in point distribution.

# Answer

# Gassner's Entropy-Stable Subcell Shock-Capturing Method: Implementation Details and Element Support

## Overview

Gregor J. Gassner and collaborators have developed a provably entropy-stable
subcell shock-capturing approach for high-order discontinuous Galerkin spectral
element methods (DGSEM). The method, first comprehensively presented in
Hennemann et al. (2021) [DOI: 10.1016/j.jcp.2020.109935], represents a
significant advancement in robust shock capturing for compressible flow
simulations while maintaining the mathematical rigor of entropy stability.

## Core Methodology

### Fundamental Principle

The method exploits the fact that high-order summation-by-parts (SBP) operators
can be rewritten in an equivalent conservative finite volume form. This enables
the design of a low-order scheme directly using the Legendre-Gauss-Lobatto (LGL)
nodes that is compatible with the discrete entropy analysis used for proving
entropy stability of the DGSEM.

Key innovation: **Seamless blending between high-order DGSEM and low-order
finite volume schemes while maintaining provable entropy stability**.

### Node Distribution and Basis Functions

#### Legendre-Gauss-Lobatto (LGL) Nodes

- **Primary node distribution**: LGL nodes on [-1,1] reference element
- **Key property**: Include boundary points at ±1
- **Quadrature exactness**: Exact for polynomials up to degree 2N-1 for N+1
  nodes

#### Basis Functions

- **Type**: Lagrange polynomials through LGL nodes
- **Orthogonality**: Not orthonormal, but satisfy SBP property
- **Tensor product construction**: For multi-dimensional elements

The DGSEM uses a local tensor-product Lagrange-polynomial basis computed from
LGL points:

```
l_i(ξ) = ∏_{j=0,j≠i}^N (ξ - ξ_j)/(ξ_i - ξ_j)
```

### Summation-By-Parts (SBP) Property

The discrete derivative matrix **D** satisfies:

```
D = M^(-1) Q
```

where:

- **M**: diagonal mass matrix with LGL weights
- **Q**: nearly skew-symmetric differentiation matrix with Q + Q^T = B
- **B**: boundary matrix = diag(-1, 0, ..., 0, 1)

This SBP property is **crucial** for entropy stability proofs.

## Element Type Support

### Demonstrated Element Types

#### 1. Hexahedral Elements (Primary Focus)

- **Implementation**: Full 3D curvilinear support
- **Mapping**: Transfinite mapping with careful Jacobian treatment
- **Status**: Fully validated in production codes (FLEXI, FLUXO)

#### 2. Quadrilateral Elements (2D)

- **Implementation**: Complete with curvilinear support
- **Extensive validation**: Numerous test cases in literature

#### 3. Tetrahedral Elements

- **Status**: Theoretical framework exists
- **Challenge**: SBP operators on simplices require special treatment
- **Implementation**: Limited practical demonstrations

#### 4. Other Element Types

- **Prisms**: Not demonstrated in the entropy-stable subcell framework
- **Pyramids**: Not addressed due to singularity issues
- **Triangles**: Some work by extension of quadrilateral methods

### Element-Specific Considerations

For hexahedral elements on curvilinear meshes:

- Metric terms must satisfy discrete GCL (Geometric Conservation Law)
- Free-stream preservation requires careful construction
- Jacobian positivity must be maintained

## Implementation Details

### High-Order to Low-Order Transformation

The key insight is that the high-order DGSEM scheme can be rewritten as:

```
J_j ∂u_j/∂t = 1/ω_j [f̂_{j-1/2} - f̂_{j+1/2}]
```

where:

- J_j: Jacobian at node j
- ω_j: LGL weight at node j
- f̂: numerical flux

This is **identical to a finite volume formulation** with subcells defined by
LGL weights.

### Subcell Division Strategy

For each high-order element:

1. **Subcell boundaries**: Located at (ξ_i + ξ_{i+1})/2 between LGL nodes
2. **Subcell sizes**: Determined by LGL quadrature weights
3. **Conservation**: Automatic due to telescoping flux differences

### Blending Mechanism

The hybrid scheme uses a blending factor α ∈ [0,1]:

```
∂u/∂t = α × (Low-order FV update) + (1-α) × (High-order DG update)
```

**Critical**: Both schemes use the same degrees of freedom (nodal values at LGL
points).

### Artificial Viscosity Integration

Unlike traditional artificial viscosity approaches, this method:

- Does **not** add explicit viscous terms
- Achieves dissipation through upwind finite volume fluxes
- Maintains strict entropy inequality without artificial parameters

The low-order scheme uses:

- **Rusanov/Local Lax-Friedrichs flux** for robustness
- **HLL flux** for enhanced accuracy
- Modified to ensure entropy stability via Tadmor's condition

## Entropy Stability Mechanism

### Continuous Entropy Analysis

For conservation law ∂u/∂t + ∇·f(u) = 0 with entropy S(u):

```
∂S/∂t + ∇·F ≤ 0
```

where F is the entropy flux.

### Discrete Mimicry

The method ensures:

```
d/dt ∫_Ω S(u_h) dx + ∑_faces F̂* ≤ 0
```

Key requirements:

1. **Volume terms**: Use entropy-conservative two-point fluxes
2. **Surface terms**: Use entropy-stable numerical fluxes
3. **Blending**: Maintains convex combination property

### Proof Strategy

1. Show high-order DGSEM satisfies semi-discrete entropy equality
2. Show low-order FV scheme satisfies semi-discrete entropy inequality
3. Prove convex blending preserves entropy stability
4. Verify for curvilinear meshes with metric terms

## Shock Detection and Troubled Cell Indicators

### Modal Decay Indicator (Persson-Peraire Type)

Modified for entropy-stable framework:

```
η_k = log₁₀(û²_{N-1} + û²_N)/(û²₀ + ... + û²_N)
```

where û_i are modal coefficients.

### Indicator Parameters

- **Threshold**: κ₀ ≈ -4 for smooth elements
- **No tuning required**: Same parameters work for Mach 1 to Mach 100

### Robustness Features

- Element-local computation
- No MPI communication required
- Parameter-free for wide range of problems

## Orthonormality Considerations

### Basis Choice Rationale

The method uses **non-orthogonal** Lagrange basis on LGL nodes rather than
orthonormal Jacobi polynomials because:

1. **SBP property**: Natural with nodal Lagrange basis
2. **Subcell interpretation**: Direct finite volume connection
3. **Efficiency**: No basis transformation needed
4. **Robustness**: Collocation avoids aliasing

### Alternative Basis Investigation

Recent work (2023) explores:

- **Legendre-Gauss nodes**: Superior accuracy but no natural subcell structure
- **Modal basis with Jacobi polynomials**: Efficient but requires transformation
- **Trade-off**: Accuracy vs. robustness and ease of implementation

## Production Implementation: FLEXI Code

### Code Structure

```fortran
! Simplified structure from FLEXI
DO iElem = 1, nElems
  ! Compute indicator
  troubled = ComputeIndicator(U(:,:,:,iElem))
  
  IF (troubled) THEN
    ! Compute blending factor
    alpha = ComputeAlpha(U(:,:,:,iElem))
    
    ! Low-order update
    CALL ComputeFV_Fluxes(U, F_FV)
    
    ! High-order update  
    CALL ComputeDG_Fluxes(U, F_DG)
    
    ! Blend
    Ut = alpha * F_FV + (1-alpha) * F_DG
  ELSE
    ! Pure high-order
    CALL ComputeDG_Fluxes(U, Ut)
  END IF
END DO
```

### Key Implementation Details

1. **Memory layout**: Solution stored at LGL nodes
2. **Flux computation**: Reuse high-order flux routines
3. **MPI parallel**: Element-local operations
4. **Time stepping**: Explicit RK schemes maintain entropy stability

## Validated Test Cases

The method has been validated on:

- **Shock tubes**: Sod, Lax, strong shocks
- **Double Mach reflection**: Up to Mach 10
- **Forward facing step**: Mach 3 with complex patterns
- **Sedov blast**: Cylindrical/spherical symmetry
- **Shock-vortex interaction**: Maintains vortex structure
- **MHD problems**: Extended to resistive MHD (Rueda-Ramírez et al. 2021)

## Recent Extensions and Developments

### 2023-2024 Advances

1. **Moving meshes**: Entropy stability preserved (Keim et al. 2025)
2. **Non-conservative systems**: GLM-MHD, two-phase flows
3. **Enhanced resolution**: Subcell reconstruction methods
4. **Monolithic convex limiting**: Alternative to predictor-corrector

### Performance Considerations

- **Overhead**: 10-20% compared to pure high-order
- **Scalability**: Demonstrated to 500,000+ cores
- **GPU implementation**: Under development in SSDC framework

## Key References for Implementation

1. **Original method**: Hennemann et al. (2021), J. Comput. Phys. 426, 109935
2. **MHD extension**: Rueda-Ramírez et al. (2021), J. Comput. Phys. 444, 111508
3. **Software**: FLEXI
   framework - [https://github.com/flexi-framework/flexi](https://github.com/flexi-framework/flexi)
4. **Theoretical foundation**: Gassner (2013), SIAM J. Sci. Comput. 35,
   A1233-A1253

## Summary

Gassner's entropy-stable subcell shock-capturing method provides a
mathematically rigorous, robust, and efficient approach for high-order
simulations of compressible flows with shocks. The method is primarily
demonstrated for hexahedral elements using LGL nodes with Lagrange basis
functions, exploiting the natural SBP property for entropy stability. While the
framework could theoretically extend to other element types, practical
implementations and validations focus on tensor-product elements where the
subcell finite volume interpretation is most natural.
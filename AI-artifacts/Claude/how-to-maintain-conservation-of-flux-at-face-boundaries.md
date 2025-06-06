# Question

What are the necessary and sufficient conditions to ensure conservation of flux
at element neighbor interfaces?

## Prompt1

Research the necessary and sufficient conditions for conservation of flux in
nodal DG finite element methods where the node points include interior and
face/edge/vertex points. Answer specifically:

1) How many face points are needed when the solution and flux polynomial are
   order P? Does this couple with the answer to (2)?
2) What technique is mandatory to conserve flux at the interface? By technique,
   I mean is some form of integral required, or can the face points simply be
   interpolated to a common surface between neighbor elements and a Riemann
   problem solution used, then corrected fluxes interpolated back to each
   neighbor's node points? Or "Lift"ed from the flux surface to the element node
   points in total from the common flux surface? What requirements would be on
   the lift operator to ensure flux conservation?
3) For conforming and non-conforming face connections between neighbor elements,
   are there different approaches required to maintain conservation of flux?
4) Are there many choices of interpolation and/or quadrature that all enforce
   conservation? What are the attributes that MUST be present in the technique
   to ensure conservation of flux at the neighbor interface?
5) What are some straightforward tests of flux conservation that can be
   performed at the element level? If there are element operators like Lift,
   Interpolation and Quadrature for face connections, what unit tests make sense
   to ensure that these operators comply with the necessary requirements for
   conservation of flux?
6) What are the performance implications of the different choices of how to
   ensure conservation of fluxes at the interfaces?

## Answer

# Flux Conservation in Nodal Discontinuous Galerkin Methods

The conservation of fluxes in nodal discontinuous Galerkin (DG) finite element
methods represents a fundamental challenge at the intersection of mathematical
rigor and computational efficiency. This comprehensive analysis examines the
necessary and sufficient conditions for achieving flux conservation across
element interfaces, drawing from seminal work by Hesthaven, Warburton, Kopriva,
Gassner, and other leading researchers in high-order numerical methods.

## Mathematical foundations for quadrature requirements

The cornerstone of flux conservation in nodal DG methods lies in the *
*fundamental Gaussian quadrature principle**: an n-point Gaussian quadrature
rule exactly integrates polynomials of degree 2n-1 or less. For polynomial
approximations of order P in DG methods, this principle dictates specific
requirements for both volume and surface integrals.

For **volume integrals**, the mass matrix entries ∫φᵢφⱼ dx involve polynomials
of degree 2P, requiring P+1 quadrature points for exact integration. However,
the stiffness matrix entries ∫φᵢ'φⱼ dx involve polynomials of degree 2P-1,
allowing P points to suffice. This distinction becomes critical when considering
computational efficiency versus accuracy trade-offs.

The **surface quadrature requirements** prove more stringent for conservation.
Face integrals of the form ∫f(u)φᵢ·n ds require at least P+1 points for
polynomial solutions of order P to guarantee flux conservation. The mathematical
proof follows from the conservation theorem requiring ∫ₓⱼ₊₁/₂ [[f(u)]] ds = 0
across interfaces, where [[·]] denotes the jump operator. Underintegration at
surfaces violates this condition even if the scheme remains stable, highlighting
that **surface quadrature determines conservation accuracy more than volume
quadrature**.

Kopriva and Gassner's spectral element analysis reveals that Gauss and
Gauss-Lobatto quadrature formulations achieve algebraic equivalence when global
polynomial interpolations are used for solutions and fluxes with appropriate
quadrature order. Their dispersion analysis shows Gauss quadrature requires
fewer points per wavelength but introduces implementation complexity through
non-diagonal mass matrices.

## Interface flux conservation techniques

The choice between direct quadrature evaluation and interpolation-based
approaches for flux conservation depends critically on problem characteristics
and desired properties. **Direct flux evaluation at quadrature points becomes
mandatory** for nonlinear conservation laws like the Euler or Navier-Stokes
equations, where flux functions are non-polynomial. Additionally, when solutions
are under-resolved relative to polynomial degree, aliasing errors can destroy
conservation and introduce non-physical oscillations.

**Interpolation-based approaches can preserve conservation** under specific
mathematical conditions. The critical requirement is that discrete operators
satisfy the summation-by-parts (SBP) property: D_SBP = M⁻¹(Q + B), where M
represents the mass matrix, Q is skew-symmetric, and B captures boundary
contributions. For conservation to hold, the identity ∑_elements ∫_Ω ∂u/∂t
dΩ + ∑_interfaces [F*] = 0 must be satisfied, requiring the jump in numerical
flux [F*] to vanish.

**Lift operators** play a crucial role in mapping surface flux differences to
volume contributions through L[F] = M⁻¹ ∫_∂K F n φ dS. Warburton established
that these operators can be precomputed for reference elements while maintaining
spectral accuracy, provided geometric transformations preserve metric
identities. However, aliasing errors in lift operators can destroy conservation
by creating spurious high-frequency modes and violating the discrete divergence
theorem.

The distinction between **quadrature-free and quadrature-based implementations**
highlights fundamental trade-offs. Quadrature-free methods offer computational
efficiency by avoiding innermost quadrature loops but only maintain exact
conservation for polynomial fluxes. For nonlinear systems, these methods
introduce flux approximation errors that compromise conservation. Conversely,
quadrature-based methods handle arbitrary flux functions and control aliasing
errors through sufficient quadrature points (typically P+2 for degree P
polynomials), guaranteeing conservation with proper quadrature order.

## Conforming versus non-conforming interface approaches

The treatment of non-conforming interfaces represents one of the most
challenging aspects of maintaining flux conservation in DG methods. **Mortar
element methods** address this challenge by introducing Lagrange multipliers on
interface boundaries, weakly enforcing both solution continuity and flux
conservation through the constraint ∫_Γ [flux] · n dΓ = 0.

For **non-conforming meshes**, mortar methods require careful quadrature
treatment. Volume and surface quadrature nodes no longer align, significantly
increasing the number of required flux evaluations. The mortar projection
operators must satisfy L²-projection properties to preserve mass locally while
handling discontinuous solution transfers. Properly constructed mortar
projections ensure the discrete conservation law d/dt ∫_Ω u dΩ + ∫_∂Ω F·n dΓ = 0
holds exactly.

**Direct interpolation methods** for non-conforming interfaces show mixed
conservation performance. They successfully conserve for linear problems with
smooth solutions where interpolation operators exactly reproduce constants and
linear functions. However, they fail for nonlinear hyperbolic problems where
interpolation introduces aliasing effects, under-resolved simulations where
gradients cannot be captured accurately, and cases with strong discontinuities
creating spurious oscillations.

The **quadrature requirements differ significantly** between interface types.
Conforming interfaces with matching polynomial orders require standard (P+1)
-point Gauss quadrature for exact flux integration. Non-conforming interfaces
with different polynomial orders must accommodate max(P₁,P₂), requiring adaptive
quadrature and careful projection between spaces. Hanging nodes from
h-adaptivity create T-junctions requiring special flux evaluation procedures and
non-uniform quadrature point distribution.

## Interpolation and quadrature rule selection

The choice between **Gauss-Legendre and Gauss-Lobatto-Legendre (GLL) quadrature
** fundamentally impacts both conservation properties and computational
efficiency. Gauss-Legendre quadrature integrates polynomials exactly up to
degree 2n-1 for n points, providing superior accuracy with fewer points per
wavelength. However, it results in a non-diagonal mass matrix requiring
expensive matrix inversion.

GLL quadrature integrates polynomials exactly only up to degree 2n-3 for n
points, falling one degree short for exact mass matrix integration—a phenomenon
known as "mass lumping." Despite this theoretical deficiency, **GLL nodes offer
critical advantages**: diagonal mass matrices enabling O(n) operations versus O(
n³) for full matrices, inclusion of boundary points simplifying boundary
condition implementation, and natural satisfaction of the SBP property for
entropy-stable schemes.

The **summation-by-parts property** proves essential for conservation and
stability. SBP operators mimic integration by parts at the discrete level
through D^T M + M D = B, ensuring discrete conservation laws mirror continuous
ones. Gassner demonstrated that Lobatto-based DG naturally possesses SBP
structure, enabling energy estimates for linear problems and entropy stability
for nonlinear conservation laws.

**Aliasing and de-aliasing requirements** become critical for nonlinear
problems. Under-integration of nonlinear terms creates spurious high-frequency
modes that can cause blow-up in under-resolved turbulent flows. De-aliasing
strategies include over-integration with higher-order quadrature (
computationally expensive), split-form methods using skew-symmetric formulations
to reduce aliasing, and modal filtering to damp high-frequency components.

## Verification procedures and testing

Rigorous verification of flux conservation requires systematic testing at
multiple levels. **Element-level conservation tests** verify that the sum of
fluxes at element interfaces equals zero for interior faces and that the DG
scheme preserves conservation at the discrete level through proper numerical
flux design.

**Unit tests for operators** must verify lift operators correctly map face
fluxes to element interiors, interpolation matrices preserve polynomial
exactness, quadrature rules integrate polynomials exactly to required degrees,
and matrix-free sum-factorization kernels match matrix-based computations.

**Free-stream preservation tests** on curved and moving meshes verify the
geometric conservation law (GCL), ensuring constant solutions are preserved
exactly on wavy and distorted grids. The method of manufactured solutions (MMS)
enables order-of-accuracy verification through analytical solutions with source
terms, using sufficiently complex solutions to exercise all governing equation
terms.

**Standard test cases** for exposing conservation errors include the
Taylor-Green vortex for energy cascade assessment, isentropic vortex for
advection accuracy on various meshes, Sod shock tube for conservation across
discontinuities, and double Mach reflection for complex shock interactions.

## Performance implications and optimization

The computational cost of ensuring flux conservation varies dramatically with
implementation approach. **Matrix-free implementations** achieve 60% of peak
arithmetic performance on modern CPUs compared to traditional matrix storage,
with sum-factorization reducing complexity from O(P^(2d)) to O(P^(d+1)) for
polynomial degree P in d dimensions. However, performance often becomes limited
by memory bandwidth (40-50% of peak) rather than arithmetic throughput.

**Memory requirements** differ substantially between approaches. Traditional DG
requires storing block diagonal matrices of size (ndof×neqs)
²×nelems—approximately 10,000 words per element for 4th-order 3D schemes.
Matrix-free approaches eliminate operator matrix storage, significantly reducing
memory footprint while maintaining accuracy through on-the-fly computation.

**Parallelization considerations** reveal that MPI ghost cell exchanges
considerably impact performance, while OpenMP demonstrates good shared-memory
scaling up to 28 cores with proper scheduling. GPU implementations require
different strategies for low versus high polynomial degrees, using specialized
kernels for low-order methods and BLAS libraries for high-order approximations.

The **trade-offs between conservation approaches** manifest in several
dimensions. Matrix-free methods provide superior computational efficiency but
may require sophisticated preconditioning. Perfect conservation may demand
additional computational overhead through flux limiting or constraint
enforcement. The optimal implementation strategy depends strongly on polynomial
degree, with tailored kernels optimal for low-order methods and BLAS-based
approaches superior for high-order approximations.

## Conclusions and best practices

Achieving robust flux conservation in nodal DG methods requires careful
attention to mathematical foundations, implementation details, and verification
procedures. **For nonlinear problems**, quadrature-based flux evaluation with
sufficient integration points remains mandatory. Surface quadrature requirements
typically dominate conservation accuracy, requiring at least P+1 points for
polynomial degree P solutions.

**For conforming meshes**, standard DG formulations with appropriate quadrature
rules naturally preserve conservation. **For non-conforming meshes**, mortar
methods provide the most robust approach, though they increase implementation
complexity and computational cost. The choice between Gauss-Legendre and
Gauss-Lobatto quadrature depends on whether accuracy or efficiency takes
priority, with GLL nodes essential for entropy-stable formulations despite their
mass lumping approximation.

**Practical implementation** should prioritize matrix-free approaches for
competitive performance, employ sum-factorization for high-order methods,
implement systematic verification through manufactured solutions and standard
test cases, and monitor discrete conservation quantities throughout simulations.
The research demonstrates that successful DG implementations require balancing
mathematical rigor with computational efficiency, using the appropriate
combination of quadrature rules, interface treatments, and implementation
strategies for the specific problem at hand.

## Prompt 2

Research the following statements, produced by Claude in answer to a prior
research query. They are inconsistent with other statements made elsewhere in
the report where mortar methods and other approaches that are interpolative are
sufficient to ensure conservation of fluxes. What does "Direct flux evaluation
at quadrature points" mean and is it actually necessary in the strict sense for
conservation of flux? Also, answer the question: Given a point distribution like
the Hesthaven Warp/Blend points that is not GLL and that doesn't have the full (
P+1)(P+2)/2 2D face point count needed for quadrature on that face, how is it
that the Hesthaven warp/blend points achieve conservation of flux for non linear
fluxes (it does) given your statements?

The choice between direct quadrature evaluation and interpolation-based
approaches for flux conservation depends critically on problem characteristics
and desired properties. *Direct flux evaluation at quadrature points becomes
mandatory* for nonlinear conservation laws like the Euler or Navier-Stokes
equations, where flux functions are non-polynomial. Additionally, when solutions
are under-resolved relative to polynomial degree, aliasing errors can destroy
conservation and introduce non-physical oscillations. *Interpolation-based
approaches can preserve conservation* under specific mathematical conditions.
The critical requirement is that discrete operators satisfy the
summation-by-parts (SBP) property: D_SBP = M⁻¹(Q + B), where M represents  
the mass matrix, Q is skew-symmetric, and B captures boundary contributions.
For conservation to hold, the identity ∑elements ∫Ω ∂u/∂t dΩ + ∑_interfaces 
[F] = 0 must be satisfied, requiring the jump in numerical flux [F] to vanish.

## Answer

# Flux conservation in discontinuous Galerkin methods: resolving apparent contradictions

The research reveals important nuances about flux conservation in discontinuous
Galerkin (DG) methods that resolve the apparent contradictions in your
questions. The key insight is that conservation and accuracy are achieved
through different mechanisms, and "direct flux evaluation at quadrature points"
provides strong guarantees but is not strictly mandatory.

## Direct flux evaluation: what it means and when it's necessary

"Direct flux evaluation at quadrature points" refers to computing the flux
function f(u) exactly at designated quadrature nodes within each element, then
using these values directly in numerical integration without intermediate
interpolation. This approach provides the strongest theoretical guarantees for
conservation, particularly for nonlinear conservation laws where polynomial
approximation of fluxes can introduce aliasing errors.

However, research demonstrates that this approach is **not strictly mandatory**
for conservation. Interpolation-based methods can preserve conservation when
they satisfy specific conditions: consistent interface treatment, sufficient
polynomial representation of the flux, and proper aliasing control. The key
requirement is that numerical fluxes at element interfaces must be computed
identically from both sides, ensuring the telescoping property that guarantees
global conservation.

For nonlinear conservation laws, the situation is more nuanced. While
interpolation can work, it faces greater challenges due to aliasing errors -
when a nonlinear flux f(u) is evaluated using degree-p polynomials, it produces
terms of degree 2p or higher, requiring more quadrature points for exact
integration. Under-integration can cause unphysical energy transfer between
different modes, potentially violating entropy conditions even while maintaining
formal conservation.

## How Hesthaven Warp/Blend points achieve conservation

The Hesthaven Warp/Blend nodal points achieve conservation through the
fundamental structure of the DG weak formulation rather than through exact
quadrature. These points are optimized for interpolation accuracy (minimizing
Lebesgue constants) rather than quadrature exactness, and they don't provide
the (P+1)(P+2)/2 face points theoretically required for exact 2D quadrature.

Conservation is maintained because of the **telescoping property** of the DG
formulation. The weak form ensures that interface fluxes cancel exactly when
summed over all elements, leaving only physical boundary contributions. As long
as numerical fluxes are evaluated consistently at interfaces (same value used
for both adjacent elements with appropriate sign), global conservation is
preserved regardless of interior quadrature accuracy.

This reveals a crucial distinction: conservation is ensured by consistent
numerical flux evaluation and the weak formulation structure, while high-order
accuracy comes from optimized interpolation points. The Warp/Blend points
provide excellent interpolation properties while the DG framework guarantees
conservation through its mathematical structure.

## Reconciling the apparent contradiction

The apparent contradiction between claims that direct flux evaluation is
mandatory versus successful interpolation-based approaches is resolved by
understanding that different methods achieve conservation through different
mechanisms:

**Direct evaluation methods** provide conservation through exact integration of
the weak form, eliminating aliasing errors by using sufficient quadrature
points. This approach offers the strongest theoretical guarantees and is most
robust for under-resolved simulations.

**Interpolation-based methods** (including mortar methods) achieve conservation
through careful construction that preserves the essential telescoping property.
Mortar methods use projection operators that maintain conservation through L²
projections, while nodal methods rely on consistent interface flux computation.

The statement that direct evaluation is "mandatory" for nonlinear conservation
laws is too strong - it's more accurate to say it provides the most robust
approach. Interpolation methods can work but require careful treatment of
aliasing errors through over-integration, filtering, or split-form formulations
that preserve kinetic energy.

## The role of summation-by-parts properties

The summation-by-parts (SBP) property D_SBP = M⁻¹(Q + B) provides a unifying
mathematical framework explaining how both approaches can achieve conservation.
Here:

- **M** is the mass matrix (quadrature weights)
- **Q** is the skew-symmetric differentiation matrix
- **B** captures boundary contributions

This property enables a discrete analogue of integration by parts, ensuring that
volume contributions to conservation cancel exactly while boundary terms
maintain proper inter-element coupling. The SBP framework shows that
conservation requires:

1. Skew-symmetry of Q (ensuring interior conservation)
2. Consistent boundary treatment through B
3. Positive-definite mass matrix M

Crucially, SBP properties can be satisfied with both direct quadrature and
interpolation approaches, provided the interpolation operators are consistent
with the underlying quadrature rule. This explains why different flux evaluation
strategies can all maintain conservation when properly constructed.

## Practical implications and recommendations

For **linear conservation laws**, interpolation-based methods work reliably with
minimal constraints. The flux functions are linear in the solution, reducing
interpolation errors and eliminating aliasing concerns.

For **nonlinear conservation laws**, the choice depends on the specific
application:

- **Well-resolved simulations**: Standard nodal DG with modest
  over-integration (1.5× minimum quadrature points)
- **Under-resolved flows**: Split-form DG methods that preserve kinetic energy
  while controlling aliasing
- **Shock-dominated problems**: Hybrid DG/finite volume approaches with robust
  limiting

The Hesthaven Warp/Blend implementation succeeds because it combines good
interpolation properties with careful numerical flux design. While it may not
achieve exact quadrature for nonlinear terms, it provides sufficient accuracy
for practical convergence while the DG framework ensures conservation.

In summary, conservation in DG methods is fundamentally about maintaining
consistency at element interfaces and preserving the telescoping structure of
the weak formulation. Direct flux evaluation provides the strongest guarantees,
but interpolation-based approaches can also maintain conservation with
appropriate design. The choice between methods involves balancing theoretical
robustness, computational efficiency, and the specific requirements of the
problem being solved.

# Prompt 3

Research the following statements about the mortar technique for ensuring
conservation of flux across non conformal element interfaces specifically for
the cases of tet/hex/prism/pyramid interfaces. Obtain the necessary
implementation detail to compose an implementation of the mortar methods in
source code in an efficient form sufficient to handle the common cases among
those elements, in particular those involving multiple triangles interfacing to
a quadrilateral for viscous and inviscid fluxes. If the implementation is
neutral WRT handling viscous/inviscid/em fluxes that is superior in that only
one method should be required to handle physical fluxes. The interest here
excludes the “how” of solving the Riemann problem, that should be a completely
separate concern from how to project fluxes to the common interface space and
back again with conservation being enforced globally and even locally if deemed
necessary.

Equal and mixed order interfaces should be included. If they differ in
implementation, call that out for a specific variant implementation for each. I
have attached the discussion to guide the research as well.

You’re absolutely right - I made an error in my explanation. Mortar methods 
are inherently asymmetric in their projection operations, and this does 
seem to contradict the telescoping principle I just described. Let me 
correct this and explain how mortar methods actually maintain conservation 
despite this asymmetry.

## The Asymmetry in Mortar Methods

You’ve correctly identified that mortar methods involve two distinct operations:

1. Forward projection: Each element projects its solution to the mortar space (
   often different for each element due to different polynomial orders or
   element sizes)
2. Backward projection: The mortar space solution is projected back to each
   element’s space (again, different for each element)
   This is fundamentally asymmetric - element A might project to a P=3 mortar
   space while element B projects to a P=4 mortar space, and the backward
   projections are similarly different.

## How Conservation is Maintained Despite Asymmetry

The key is that mortar methods maintain conservation through a different
mechanism than simple telescoping. They ensure conservation by satisfying a
discrete compatibility condition:

∫_Γ [F*]_mortar · v_A dΓ = -∫_Γ [F*]_mortar · v_B dΓ

for all test functions v_A and v_B in the respective element spaces. This is
achieved through:

1. L² projections that preserve integral quantities: The projection operators
   are designed so that integral conservation is maintained even though
   pointwise values differ
2. Careful construction of the mortar flux: The numerical flux F* is computed in
   the mortar space using a single consistent evaluation, then projected back to
   each element space in a way that preserves the integral of the flux
3. Compatibility of quadrature rules: The quadrature rules on each side must be
   compatible with the mortar space quadrature in a specific sense

## The Subtle Distinction

The conservation property in mortar methods comes from preserving the integral
of the flux across the interface, not from having identical pointwise values.
This is a weaker but sufficient condition for conservation.
In mathematical terms, if P_A and P_B are the projection operators to elements A
and B from the mortar space:

* P_A(F) ≠ P_B(F) (pointwise)
* But ∫ P_A(F) · v_A = -∫ P_B(F) · v_B (in the weak sense)

This research requires understanding the state of practice for mortar
interfaces. If common and necessary practices form a mortar space a certain way,
say using a quadrature surface of a specific degree in order to satisfy the
conservation or accuracy (or stability) requirement, this should be discovered
as part of the research.

The resulting report should contain enough detail for one, maybe two
implementations of mortar interfaces. Possibly one for equal order polynomial
neighbors and another for mixed order. It should be opinionated WRT an
implementation approach that mirrors the most successful approaches used for
mixed element production DG FR methods for CFD and MHD. The target application
is to a method library that will provide element bases for the four 3D element
shapes, and the mortar interface will need to support all of them.

Another detail to incorporate: the element definitions have not yet standardized
on a node point choice. I can choose from any of the most successful types for
CFD flows with strong shocks and iLES, including GLL/LGL points, Hesthaven
Warp/Blend or others as needed. If it is discovered as part of the research that
the volumetric node choice is important to the success of the flux continuity
scheme, I am open to that being part of the implementation AS AN ASSUMPTION,
meaning make the overt statement that the implementation should assume the use
of one of the volume node point choices, and all of the implementation will
proceed with this assumption. This becomes important as to the origin of the
flux (solution node point directly or interpolated), etc.
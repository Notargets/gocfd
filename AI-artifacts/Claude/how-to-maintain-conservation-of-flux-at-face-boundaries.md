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
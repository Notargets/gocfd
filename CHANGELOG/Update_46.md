# After Report: No profound new technique for flux projection is advisable

After reviewing the strategy outlined in Update 45, I've concluded that there
aren't any new benefits that would arise from a new technique to project the 
vector flux onto the RT edge nodes after all. Instead, I've begun revisiting
various other aspects of the methods used right now and the test cases 
implemented to validate them.

## Specifically, why not use a new flux projection technique?

The primary approach of the Galerkin Finite Element method is to capture a 
solution field as a polynomial function on the defined element. In this case,
we are using the finite element to represent the conserved field values like 
mass, momentum and energy. The divergence of the flux of these variables is 
a calculated quantity from those fields, the calculation of which involves 
derivatives that degrade the polynomial order by one degree.

The Raviart Thomas element is used to fix the reduction of polynomial 
order by carrying the computed flux of the scalar fields at one polynomial
degree higher than the scalar element, then the derivatives needed to compute
divergence can be performed at the order of the scalar element.

The step involving transfer of the vector flux to the points that define the 
RT element should be done using the definition of the polynomial for the 
scalar fields. Interpolation is the right approach for this step, as it 
correctly represents what we're capturing, a scalar polynomial field sampled 
at the nodes defining the RT element.

## Are there alternatives to interpolation? What about projection?

While projection can be used, there are many challenges associated with it 
that muddy the nature of the resulting method. In particular, projection 
uses integration of the scalar field over the element, then a redistribution 
of the values of the polynomial over the original nodes that define the 
polynomial. This is an ill-defined problem requiring techniques like least 
squares fitting to re-distribute values across the element. Again, it seems 
incongruous to have to "smooth" the element's polynomial representation by 
projection arbitrarily.

## So, why is the method so unstable at almost all orders right now?

After reviewing the existing tests, which do cover many accuracy tests of 
the divergence and derivative calculations for polynomial fields, it's clear 
that the current implementations are, in-fact, accurate. Shock tube results 
and prior isentropic vortex calculations are very accurate. The principle 
problem we have is stability when fields aren't smooth.

One thing that has grabbed my attention is the formulation of the Raviart 
Thomas element's basis functions. While the scalar element is implemented 
using orthogonal polynomials along the lines of Hesthaven and Warburton, the 
RT elements have been defined using base cartesian polynomials. It is 
possible that the disconnect between these two field representations is 
causing amplification of non-smooth waves as interactions between the scalar 
and RT elements.

## Current Plan: New RT element formulation to match the scalar elements

I'm planning to implement the orthogonal basis functions identical to the 
scalar element for the interior nodes of the RT element. The edge functions 
of the RT element will also use those functions, multiplied by the base 
normal vectors on the edges, with special care to ensure the polynomials are 
smooth along edge traversals.
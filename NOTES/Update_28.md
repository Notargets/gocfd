### Update: July 7 2020):

Researching the use of the Raviart-Thomas finite element to represent the numerical flux and divergence.

I've found a lot of good material on these elements and the Raviart-Thomas element, but I'm struggling to find a simple connection between the convenient techniques used for Lagrange elements as expressed in the Hesthaven text and non-Lagrangian elements. In particular, the Lagrange elements have the property that the interpolation operator and the modal representation are connected simply through the Vandermonde matrix. This enables the easy transformation between function values at nodes within the element and interpolations and derivatives. If we could use a Raviart-Thomas element similarly to a Lagrange element, we'd just compose a Vandermonde matrix and proceed as normal to compute the divergence, but, alas, I've struggled to find a way to do this.

I have found the open source element library "Fiat", part of FEniCS, that implements RT elements and will produce a Vandermonde matrix for them, but I'm not convinced it's useful in a way similar to Lagrange elements. Also, I found that the ordering of the polynomial expansions are strange in the RT elements in Fiat - not sure it matters, but it's not a straightforward "smallest to largest" ordering, which brings into question how to use that Vandermonde matrix in general.

So - I'm back to basic research on the meaning of finite elements and how to properly represent the flux on Raviart-Thomas elements. Papers I'm reading can be found [here](research/convergence_and_fluxes).


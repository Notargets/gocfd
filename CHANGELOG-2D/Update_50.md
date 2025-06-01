# Redefinition of the Raviart Thomas Element Construction (Part 4)

## Fixed Polynomial basis - Now using orthogonal 2D Lagrange polynomials
The interior basis functions now use the Jacobi 2D basis but are implemented 
using the orthogonal Lagrange polynomials. This is done by forming the 
Vandermonde matrix inverse to obtain the Lagrange interpolating coefficients 
for each member of the basis. This works perfectly to define orthogonal 
polynomials for each of the j-th interior basis functions.

## Fixed Ervin basis - edge basis functions now used unit edge normals
The problem with the Ervin basis using non-normal, non-fixed edge 
polynomials is that the interior basis contributes tangential flux to the 
edge basis, which is not appropriate for the RT element. Now the edge basis 
functions use the unit normal vectors multiplied by an order P 1D polynomial.

Note that the interior basis functions for the Ervin basis used a fixed 
vector in each of two types E4 and E5 that, when a dot product is formed 
with each edge normal, is zero. That ensures that there is no interior 
contribution to the normal flux on the element edge.

## Interpolation of the Flux vector field onto the RT element
Since we know the values of the Flux vector field at every basis DOF for the RT
element, we can use interpolation to derive the constants that represent the
field on the polynomial basis. This is described in a ChatGPT session
[The_Interpolation_Equation_For_RT.pdf](../ChatGPT/The_Interpolation_Equation_For_RT.pdf)

## Flux projection onto the RT element
To *project* the flux onto the RT element, we need to solve an integral equation
that is described via a set of prompts to ChatGPT and a nice markup document
I've put together in
[The Projection Equation.pdf](../ChatGPT/The_Projection_Equation.pdf).

To do a projection, we would devise a set of quadrature weights for the RT 
element to enable solving the integration needed for projection. Then we can 
construct the basis evaluation equation using the quadrature weights, and 
then finally solve the equation to determine the coefficients that map the 
vector flux onto the RT element basis. With the "mass matrix" from this 
equation, we can construct the divergence matrix that will enable a simple 
way to calculate divergence using the flux vector from the Lagrange element.

Projection is only needed if interpolation does not suffice.
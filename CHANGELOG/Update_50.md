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

## Researched the proper construction of the Flux projection onto the RT element
It turns out that in order to properly construct the Divergence relationship 
with the RT element basis, we need to solve an integral equation that 
projects the vector flux function onto the element basis. This is described 
now via a set of prompts to ChatGPT and a nice markup document I've put 
together in
[The Projection Equation.pdf](../ChatGPT/The_Projection_Equation.pdf).

The next steps are to put together a set of quadrature weights for the RT 
element to enable solving the integration needed for projection. Then we can 
construct the basis evaluation equation using the quadrature weights, and 
then finally solve the equation to determine the coefficients that map the 
vector flux onto the RT element basis. With the "mass matrix" from this 
equation, we can construct the divergence matrix that will enable a simple 
way to calculate divergence using the flux vector from the Lagrange element.
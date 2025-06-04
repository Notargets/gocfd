# Prompts for the Claude design specifications

## 3D Element construction

### Primary prompt
I need a detailed implementation plan for finite element polynomial bases for 
each of the 3D element types tet (PKD), hex (sum factorization using 
Persson lines), prism (sum factorized in 2 directions), pyramid (yuck) for 
use in CFD / MHD simulations of LES with shocks and other contact 
discontinuities. I will be solving the physics equations in differential 
form and will need very efficient, accurate and stable difference operators 
for computing divergence and gradient in nodal form, and will also need to 
access modal for filtering and limiting.

Write a very detailed development plan for each element type that specifies 
each set of code to be written and the unit tests (golang) that will verify 
each step. Write the plan in such a way that the development is done 
incrementally. Functions for doing matrix operations, 1D and 2D Jacobi 
polynomials and Vandermonde matrices are provided. For instance, the PKD 
basis for tets should be broken down into steps like: 1) implement 
koornwinder basis functions for r,st. Unit test  for interpolation accuracy. 
Unit test for orthonormality of basis. 2) implement derivative matrices for 
(1). Unit test for polynomial order P monomial derivative accuracy with sum 
of monomial terms in each direction. 3) implement affine transform of 
coordinates for all element types. Unit test X,Y,Z based polynomials for 
stretched/distorted element.

The instruction should be written such that claude.ai opus 4 should be able 
to develop each element basis using the intructions in golang code in 
increments, each increment checked by the user before the next increment is 
initiated.

#### Supplements

##### 1
Research fortunato and persson papers on line based basis for tensor 
product elements and revisit the plan artifact to properly reflect how to 
implement that basis for hex and prism elements

##### 2
Insert the persson paper citations into the plan artifact.

Research the pkd basis definitions to find the specific mass matrix based 
orthonormalization formulas and dubiner normalization formulae. Update the 
plan an place the citations there.

Update the plan to include implementing the mode ordering slice for each 
element type so that modes of a given order can be accessed conveniently.

Remove from the plan extraneous info about cfd and mhd development and time 
schedules. Leave only information about specific element development tasks 
with details and specifics on what to implement and tests that prove each 
step is correct.

##### 3
Include in the plan williams-shunn-Jameson symmetric quadrature rule points 
for order up to 7 for tetrahedra with citations included in the plan.

Locate optimum symmetric quadrature rule points for the other element types 
that have the same features as the wsj points, namely the colocation of 
points with the node points for an order P element such that the quadrature 
will exactly integrate order 2 P.

## Face quadrature, Lift operator construction, non-conforming interface mortar

### Primary prompt

Write a section of documentation targeted toward describing the 
implementation detail of face quadrature points and lift operator 
construction for a 3D element library written in golang. Do not write code.
Write in document form sufficient for Claude.ai Opus 4 to be able to 
implement the code in an incremental fashion, as it will be added to an 
existing set of code that implements the various element geometry bases for 
tet (PKD), hex (line based basis), prism (hybrid line and collapsed), pyramid.

Include citations  as needed to determine the formulae necessary to 
implement the lift operators that are conservative in nature. Break the 
description into one section per each element type. Describe the work in 
incremental steps, each one with accompanying unit tests that prove that 
step is correctly implemented. For each element type, describe the steps to 
construct the mortar for each face type. If there is a shared mortar type 
used to do the calculation of the Riemann flux used by all element faces, 
describe that in an addition section with incremental instructions on 
building it with unit tests for each increment.
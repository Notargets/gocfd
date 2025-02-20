# Redefinition of the Raviart Thomas Element Construction (Final)

## The reworked RT element for 2D simplex triangles is complete

The new element completes testing through order RT7 and demonstrates 
quadratic convergence of accuracy for divergence in P for Sin/Cos vector fields.

The new element uses Jacobi polynomials in Lagrange basis form for the 2D 
and the 1D polynomials in the basis.

## Effect on Gibbs issues in CFD solver: none

Replacing the existing RT element with this one had no effect on the issues 
with the CFD solver. This points to some other issue in the method 
formulation that is causing the Gibbs problems with shock waves, along with 
the known issues with polynomial orders greater than 2, and even order 1.

At this point, I suspect the complexity of the CFD solver implementation has 
resulted in some bugs or issues that I can't find easily.

## Next step: Rebuild a new CFD solver from scratch - much better and simpler

This time, rather than using complex edge structures designed to minimize 
memory usage and parallel processing, I'll do something much much simpler to 
work with. Instead of dedicated edge data structures, I'll work to make the 
edge manipulations centered on the elements.

## Implicit solver is needed for complex / viscous solutions

We need to be able to use Implicit time advancement in order to deal with 
highly skewed meshes, which are always needed when we have high Reynolds 
number viscous flows. Mesh sizes near walls in these flows are often 1E6 
times as small as those used in flow away from walls in order to properly 
model the boundary layer, heat transfer and other phenomena. This causes the 
time steps needed to get too small to make for practical solutions. An 
implicit solver fixes this issue in large part, in exchange for using a lot 
more memory and compute time per iteration.

## Need: A sparse block matrix factorization / solver

There isn't a current block matrix (where each matrix position is itself a 
dense matrix) library that implements sparse matrices and also implements a 
solver / factorization method set like GMRES or LU.

I've found an existing sparse block matrix library that implements large 
scale matrix->vector multiplication, so I'm going to implement GMRES using 
it. This will allow for solution of large sparse block systems efficiently 
in Golang, which is apparently missing RN.
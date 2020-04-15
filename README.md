# gocfd
Awesome CFD solver written in Go 

## An implementation of the Discontinuous Galerkin Method for solving systems of equations

### Example: Maxwell's Equations solved in a 1D Cavity
![](images/Maxwell1D-cavity0.PNG)

This is one of the two model problem solvers currently implemented, a 1D Maxwell's equations solution in a metal cavity with a change of material half way through the domain. The initial condition is a sine wave for the E (electric) field in the left half of the domain, and zero for E and H everywhere else.

Unlike the advection equation model problem, this solver does have unstable points in the space of K (element count) and N (polynomial degree). So far, it appears that the polynomial degree must be >= 5 for stability, otherwise aliasing occurs, where even/odd modes are excited among grid points.

In the example pictured, there are 80 elements (K=80) and the polynomial degree is 6 (N=6).
![](images/Maxwell1D-cavity.PNG)
![](images/Maxwell1D-cavity2.PNG)

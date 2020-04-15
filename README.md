# gocfd
Awesome CFD solver written in Go 

## An implementation of the Discontinuous Galerkin Method for solving systems of equations

### Model Problem Example #1: Advection Equation
![](images/Advect1D-0.PNG)

The first model problem is 1D Advection with a left boundary driven sine wave.

In the example pictured, there are 80 elements (K=80) and the element polynomial degree is 5 (N=5).
![](images/Advect1D-1.PNG)

### Model Problem Example #2: Maxwell's Equations solved in a 1D Cavity
![](images/Maxwell1D-cavity0.PNG)

Another model problem is a 1D Maxwell's equations solution in a metal cavity with a change of material half way through the domain. The initial condition is a sine wave for the E (electric) field in the left half of the domain, and zero for E and H everywhere else. The E field is zero on the boundary (face flux out = face flux in) and the H field passes through unchanged (face flux zero), corresponding to a metallic boundary.

Unlike the advection equation model problem, this solver does have unstable points in the space of K (element count) and N (polynomial degree). So far, it appears that the polynomial degree must be >= 5 for stability, otherwise aliasing occurs, where even/odd modes are excited among grid points.

In the example pictured, there are 80 elements (K=80) and the element polynomial degree is 5 (N=5).
![](images/Maxwell1D-cavity.PNG)
![](images/Maxwell1D-cavity2.PNG)

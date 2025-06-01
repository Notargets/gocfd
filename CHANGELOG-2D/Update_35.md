## Update: [10/12/21]
I've implemented two different gradient calculations to compute the dissipation and verified they function correctly.

The dissipation is not able to suppress shock-induced (Gibbs) oscillations as implemented, and I think it's clear that this is
due to an effect described by Persson involving the continuity of the gradient operator. In the picture below (2nd order, Roe
flux, M=0.8), you can see the base of the shock wave above the NACA0012 airfoil. As the shock strengthens, the X derivative of
X momentum is growing and is clearly discontinuous between elements, despite being computed from a C0 continuous X momentum
field (accomplished using the RT element). For the sub-element shock capturing to work properly, Persson suggests that the
gradient must be smooth across elements; otherwise, the jumps between elements are strengthened instead of being diminished as
needed. The solution described is to use the "LDG" method for computing 2nd order derivatives as described by [Cockburn and
Shu](../research/filters_and_flux_limiters/cockburn-shu-LDG-second-order-terms.pdf) in 1999.

Next, I plan to implement 2nd order derivative continuity along the lines of the LDG method, which will also support the next
steps in implementation of viscous equations (Navier Stokes). Hopefully, we'll kill two birds with one stone: sharp shock
capturing at high order accuracy and viscous solutions!

![](../images/discontinuous-gradient-in-shock.PNG)


[Back to Index](../CHANGELOG-2D.md)

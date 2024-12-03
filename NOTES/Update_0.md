# gocfd
Awesome CFD solver written in Go

| NACA 0012 Airfoil at M=0.3, Alpha=6, Roe flux, Local Time Stepping | M=0.5, Alpha=0, Roe Flux, 1482 O(2) Elements, Converged |
|:------------------------------------------------------------------:|--------------------------------------------------------:|
|               ![](images/naca12_2d_m0.3_a6_roe.gif)                |                 ![](images/naca12_2d_m0.5_aoa0_Roe.PNG) |

|                           Density                            |                            X Momentum                             |                  Density                   |
|:------------------------------------------------------------:|:-----------------------------------------------------------------:|:------------------------------------------:|
| ![](images/render-mesh-isentropic-vortex-initial-zoom-7.PNG) | ![](images/render-mesh-isentropic-vortex-initial-zoom-7-rhoU.png) | ![](images/vortex-1-2-4-7-lax-cropped.gif) |

====
## Update: [11/2/24]

After a three year hiatus I'm back!

I had previously confirmed that the oscillations at orders greater than P=0 is due to interpolation through
discontinuities like shock waves within elements. I have implemented artificial dissipation and a limiter, but have not
yet been able to get a workable system with monotone solutions that converge predictably.

Enter a new test case that eliminates boundary conditions and other complexities: a 15 degree wedge with an incoming
shock wave. I ran the case with the existing options to control oscillations and show solutions below for P=0, P=2 and
P=4 for this case at M=2 with 948 triangle elements. 

The Barth/Jespersen limiter is able to control the oscillations, smearing the shock over 3 or so triangle elements, but
the solutions don't converge properly due to the widely reported "buzzing" of residual as the shock switch oscillates.

The Persson C0 dissipation operator is better behaved and the solutions converge well. THe shock is resolved within the
triangle elements instead of smearing over multiple elements. The issue with Persson is that the parameter Kappa varies
by case and the runtime penalty is very large, as we have to solve for a higher order field to get the dissipation
values, similar in complexity to implementing a viscous flow field.

|                    No Limiter<br/>P=0                    |                            P=2                             |                                P=4                                 |
|:----------------------------------------------------------:|:----------------------------------------------------------:|:------------------------------------------------------------------:|
|   ![](images/M=2-15deg-wedge-P=0-converged-nofilter.PNG)   |   ![](images/M=2-15deg-wedge-P=2-converged-nofilter.PNG)   |                        No Solution Possible                        |
|                Barth Jespersen Limiter<br/>                |                                                            |                                                                    |
|   ![](images/M=2-15deg-wedge-P=0-converged-nofilter.PNG)   |   ![](images/M=2-15deg-wedge-P=2-BarthJ-K5-buzzing.PNG)    |       ![](images/M=2-15deg-wedge-P=4-BarthJ-K5-buzzing.PNG)        |
|           Persson C0 Continuous Dissipation<br/>           |                                                            |                                                                    |
| ![](images/M=2-15deg-wedge-P=0-converged-perssonC0-K5.PNG) | ![](images/M=2-15deg-wedge-P=2-converged-perssonC0-K5.PNG) | ![](images/M=2-15deg-wedge-P=4-converged-perssonC0-K5.PNG) |

## Update: [11/28/21]

I've now continued testing and transonic cases are also working for P!=1. What's interesting now is that for P=0,
the transonic solutions and shock tube solutions are wiggle free / monotone. Below is a transonic airfoil solution at
Mach = 1.0 and 0.8, converged, on a mesh with 42,500 elements. The shock wave is finely resolved and is monotone. The convergence
was also far more monotone and rapid than the higher order solutions.

The biggest difference between the 0th order solutions and higher order is the interpolation of the flux from the interior
nodes. In the P=0 case, there is only one interior point and the flux is simply copied from the interior to the edge, as
is typical of second order finite volume schemes. Reconstruction of the shared flux value at the edge is then handled via
solving the 1D normal Riemann problem at the edge at first order. The resulting scheme should converge at 1st order, or
P+1.

I think at this point, the wiggle problem I'm experiencing is due to the non-monotonicity of the interpolation of the flux.
New extrema are being created, and my instinct is to eliminate the extrema by limiting the interpolated values to the min
and max of the values in the interior. This would remove "extrapolation artifacts", but what is the impact on numerical
accuracy?

| P=0, NACA 0012, Mach = 1.0, Alpha = 2, Interpolated Flux, Converged |
|:-------------------------------------------------------------------:|
|            ![](images/interpolateFluxNotQ-M=1.0-P=0.PNG)            |

| P=2, NACA 0012, Mach = 0.5, Alpha = 2, Interpolated Flux, Converged |     P=0, Mach = 0.8, Fine Mesh, Converged     |
|:-------------------------------------------------------------------:|:---------------------------------------------:|
|              ![](images/interpolateFluxNotQ-M=0.5.PNG)              | ![](images/interpolateFluxNotQ-M=0.8-P=0.PNG) |


## Discontinuous Galerkin Method for solving systems of equations - CFD, CEM, ... hydrodynamics-fusion (simulate the Sun), etc! 

##### Credits to:
- Jan S. Hesthaven and Tim Warburton for their excellent text "Nodal Discontinuous Galerkin Methods" (2007)
- J. Romero, K. Asthana and Antony Jameson for "A Simplified Formulation of the Flux Reconstruction Method" (2015) for the DFR approach with Raviart-Thomas elements

### Objectives

1) Implement a complete 3D solver for unstructured CFD (and possibly MHD) using the Discontinuous Galerkin (DG) method
2) Optimize for GPUs and groups of GPUs, taking advantage of the nature of the natural parallelism of DG methods
3) Prove the accuracy of the CFD solver for predicting flows with turbulence, shear flows and strong temperature gradients
4) Make the solver available for use as an open source tool

It is important to me that the code implementing the solver be as simple as possible so that it can be further developed and extended. There are other projects that have achieved some of the above, most notably the [HiFiLES](https://hifiles.stanford.edu/) project, which has demonstrated high accuracy for turbulence problems and some transonic flows with shock waves and is open source. I personally find that C++ code is very difficult to understand due to the heavy usage of indirection and abstraction, which makes an already complex subject unnecessarily more difficult. I feel that the Go language makes it easier to develop straightforward, more easily understandable code paths, while providing similar if not equivalent optimality and higher development efficiency than C++.  

### Why do this work?

I studied CFD in graduate school in 1987 and worked for Northrop for 10 years building and using CFD methods to design and debug airplanes and propulsion systems. During my time applying CFD, I had some great success and some notable failures in getting useful results from the CFD analysis. The most common theme in the failures: flows with thermal gradients, shear flows and vortices were handled very poorly by all known usable Finite Volume methods.

Then, last year (2019), I noticed there were some amazing looking results appearing on Youtube and elsewhere showing well resolved turbulent eddies and shear flows using this new "Discontinuous Galerkin Finite Elements" method...

##
Update: [12/01/21]

I've implemented a min/max limiter that restricts values interpolated from the
interior solution points to the edges in a couple of different variations. When the limiter is active, no interpolated
edge values are higher or lower than the interior, All values used derive from the solution of the fluid equations in
the interior nodes. This eliminates extrapolation overshoot at edges known as the Runge phenomenon. What I observed is
that the Gibbs instability became more pronounced, with the "wiggles" appearing throughout the interior at higher
frequencies, which confirms that the overall problem is actually in the interior polynomials, not the interpolation to
the edges, or possibly a coupled version of the edge and interior.

I did find a description of a very similar project using the Flux Reconstruction method of Huynh [Vandenhoeck and Lani, "Implicit High-Order Flux Reconstruction Solver for High-Speed Compressible
Flows"](research/flux_reconstruction/vandenhoeck-coolfluid-solver.pdf), which is a more complicated process with an
identical resulting numerical scheme (on quadrilaterals/tenso product bases). Vandenhoeck notes that the Gibbs
instability increases when the artificial viscoscity length scale is "too small", and in their case he reported that
the length scale should be large enough to spread the shock over "a couple of cells" so as to avoid the instabilities
inside the cell consequent from the intra-cell capture of the discontinuity.

I think the next step I'll take is to modify the artificial viscoscity to parameterize it similarly to Vandenhoeck, then
run some experiments to validate stability when capturing shock waves over a couple of elements.

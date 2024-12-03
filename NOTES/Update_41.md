## Update: [11/27/21]

I've implemented interpolation of fluxes to the edges, replacing the interpolation of solution values followed by
computation of the flux. It's currently an option for Lax and Roe fluxes. Below is a 2nd order converged solution using
the new approach, showing very smooth solution contours.

It's definitely a more computationally expensive route, but it seems to be "correct", in that there are fewer interpolation
errors evident in the solution. The results are substantially different in that we get:
1) stable and convergent solutions at subsonic Mach for orders 0,2,4, unstable at P=1, likely due to interpolation overshoot
2) smoother and more realistic solution contours without the edge defects

I think it's likely this is a better formulation for the edge computations, but we need something that will eliminate
the spurious new minima/maxima being introduced by the flux interpolation. TVD/ENO concepts come to mind, where we
limit the interpolated values so as not to introduce new minima/maxima, but I think it's important to formulate such an
operation so that it doesn't create aphysical effects, especially in the time accurate solver.


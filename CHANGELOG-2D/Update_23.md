## Update: (Nov 24 2020):

The 2D Euler solver now works! Yay!

I've tested the Isentropic Vortex case successfully using a Lax flux and a Riemann BC backed by an analytic solution at the
boundaries. Unlike the DG solver based on integration in Westhaven, et al, this solver is not stable when using a Dirichlet
(fixed) BC set directly to the analytic solution, so the Riemann solution suppresses unstable waves at the boundary. I
haven't calculated error yet, but it's clear that as we increase the solver order, the solution gets much more resolved with
lower undershoots, so it appears as though a convergence study will show super convergence, as expected.

The solver is stable with CFL = 1 using the RK3 SSP time advancement scheme. I'll plan to do a formal stability analysis
later, but all looks good! The movie below is 1st Order (top left), 2nd Order (top right), 4th Order (bottom left), and 7th
Order (bottom right) solutions of the isentropic vortex case with the same 256 element mesh using a Lax flux. We can see the
resolution and dispersion improve as the order increases.

![Simulation Visualization](../images/vortex-1-2-4-7-lax-cropped.gif)

You can recreate the above with ```gocfd 2D -g --gridFile DG2D/vortexA04.neu -n 1```, change the order with the "-n" option.
"gocfd help 2D" for all options.

[Back to Index](../CHANGELOG-2D.md)

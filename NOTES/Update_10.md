### Update: [9/18/21]
Update: The data are in - the new Roe-ER flux *is* faster to compute and has good characteristics, so it's a good thing  and will be useful for turbulence capturing and strong shock applications later. However, tests showed that using either *Roe* or *Roe-ER* flux calculations crashed due to odd-even instability with shocks past maybe Mach 0.6 on the airfoil test. I'm thinking that we need to stabilize the flux at the border of the element using higher order approximations and clipping them using MUSCL/TVD approach. This is similar to a typical J->J+1 structured approach but using "wiggle simulation" for the flux derivatives at the boundary. I wonder also if I'll need to increase the derivative degree at the boundary along with the overall order of the element? It's very do-able, though I think I'd want to do a more analytic derivation for sampling higher order fields within elements.

Removing the "wiggles", part X: There are still odd-even instability modes being triggered by shocks, which for now I'm assuming are being introduced by the flux transfers between elements. It's possible that the discontinuities from shocks are landing in the polynomial basis and amplifying there, but first I'm going to eliminate the flux transfer amplifications before doing anything inside the element.

I'm planning now to implement a TVD scheme at 2nd order for the flux transfer as follows:
1) Interpolate the Q field to the edge using the element polynomial basis (same as now)
2) Interpolate an additional value close to the edge, say 0.01 of the edge length close
3) Use the two values on either side of each edge to construct a MUSCL with TVD to obtain the edge flux

This approach should introduce damping to oscillatory modes crossing the boundary, while slightly improving the accuracy of the interpolated flux.


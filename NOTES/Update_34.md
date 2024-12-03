## Update: [10/3/21]

I've implemented Laplacian artificial dissipation that tracks shock-induced instabilities using the Lagrangian solution
element to compute the flux derivatives and also for the divergence of the dissipation field. The method works well for 1st
order calculations, with sharp shock resolution and fast convergence, though the intra-element shock resolution is marked with
a significant discontinuity with the edges of the shock capturing cell. When the shock aligns with the edge, the result is a
near-perfect shock capture, but when the shock is not on the edge, the intra-cell solution has spurious internal oscillations.

Based on these results, I'm implementing what seems a simpler and more appropriate method: The dissipation flux will be
calculated for all points within the Raviart-Thomas element used for the calculation of the fluid flux divergence. The
resulting dissipation flux [epsR, epsS] is then added to the physical fluid flux [fluxR, fluxS] and then the divergence of
the combined flux is calculated using the RT element divergence operator. I'm anticipating that this approach should provide
added stability to the flux field, which is an (N+1)th order polynomial field that overlays the (N)th order solution field. In
the prior approach, the (N)th order divergence of the artificial dissipation is added to the (N+1)th order derived physical
divergence, and I think that is creating aliasing errors that are unstable.

[Back to Index](../NOTES_Index.md)

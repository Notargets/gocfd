## Update: [12/01/21]

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


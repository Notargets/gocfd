## Update: [11/15/21]

New theory: it's the interpolation of the solution values to the edges that is the culprit - instead,
I'll interpolate the flux from the solution points to the edge along with the solution values, then compose the Roe
flux as an average of the interpolated L/R fluxes plus terms arising from the interpolated solution values for the
Riemann problem at the edge.


## Update: Nov 18 2020):

The 2D Euler solver is functionally complete now, although lacking (many) boundary conditions, the Roe and Lax flux
calculations and testing/validation. I want to display the solution graphically to debug.

Next up: graphics to display the solution on the full set of polynomial nodes. Graphics will be key to debugging and
answering questions like: "is it better to interpolate flux or solution values, then compute flux?" My plan is to export
the triangulation of the RT or Lagrange elements to the graphics plotting in AVS. One question to address: how to handle
the corners of the triangles? Currently, the corners are not present in the RT or Lagrange elements - we only have edges up
to the corner and interior points. I think that I would rather average the two closest edge points to create the corner
value than interpolate corners - my rationale is that I don't want to create artificial extrema in the solution while
plotting, it can distract from the actual, and likely problematic extrema I'm seeking to plot.


[Back to Index](../CHANGELOG-2D.md)

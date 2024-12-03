### Update: [8/28/21]

Working on adding routines to invert and solve block matrix systems so that I can implement an implicit time advancement scheme.

Implicit fluid dynamics schemes all use the flux jacobian of the underlying equations, which yields a system of equations that has "blocks" consisting of the 4x4 (2D) or 5x5 (3D) matrices representing each nodal point. This requires that we have a way to invert the system matrix, or use a solver approach to get the final update in the time advancement scheme.

The work I'm doing now will add a capability to efficiently store the full block system matrix and use it to solve for the time update.


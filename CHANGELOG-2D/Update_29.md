## Update: (Aug 18, 2021):

Investigating algorithmic methods to speed up convergence to a steady state solution. These techniques can also be applied
to unsteady/time varying solutions by solving incremental time advancement sub-problems.

There are many algorithms in the literature, including:

#### Node implicit
Implicit time advancement within each element, which should improve the convergence rate as the order of the elements
increases.

#### Line Gauss-Seidel Implicit
Improves global implicitness by solving advancing lines (2D) or planes (3D) through the field. This technique has been shown
to greatly improve convergence of viscous problems with cells packed tightly to surface geometry.

#### Multigrid
Propagates low frequency changes rapidly through the finest mesh using nested coarse meshes.

#### Preconditioners
Used to remove stiffness in propagating high-frequency changes where the difference in wave speeds is large. Examples include
very low speed flows where the acoustic wave speed and sonic wave speeds differ greatly, or viscous problems where the
viscous eigenvalues are very small compared to sonic waves.

[Back to Index](../CHANGELOG-2D.md)

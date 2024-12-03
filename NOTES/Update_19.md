## Update: Nov 5 2020):

Progress on the 2D Euler equations solution! There are now unit tests showing that we get zero divergence for the
freestream initialized solution on a test mesh. This includes the shared face normals and boundaries along with the
solution interpolation. The remaining work to complete the solver includes boundary conditions and the Roe/Lax Riemann
flux calculation at shared faces - each of which are pure local calculations.

I'm very happy with the simplicity of the resulting algorithm. There are two matrix multiplications (across all elements),
a matrix multiplication for the edge interpolation and a single calculation per edge for the Riemann fluxes. I think the
code is easily understood and it should be simple to implement in GPU and other parallel systems.

[Back to Index](../NOTES_Index.md)

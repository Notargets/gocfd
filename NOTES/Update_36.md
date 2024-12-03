### Update: (May 9, 2020): Direct Flux Reconstruction implemented for 1D Advection, moving to implement for the other two 1D model problems

DFR works for Advection (-model 3) and seems to improve accuracy and physicality.

The primary differences between using DFR and traditional Nodal Discontinuous Galerkin:
1) Instead of using only the primitive variables in the right hand side, we use the Flux directly for a hyperbolic problem
2) The flux is a globally composed variable and is differentiated after being "reconstructed" to be C(0)
3) The reconstruction technique follows Jameson(2014) in that we coerce the flux values of each face to be consistent values at each face, then use the same Np Gauss-Lobato nodes and metrics to develop the derivative(s) of flux

For (3) in this case, I used the same Lax-Friedrichs flux calculation from the text, then averaged the values from each side of each face together to make a single consistent face value shared by neighboring elements.


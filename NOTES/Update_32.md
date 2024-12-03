### Update: June 2, 2020): 

Implemented a smooth solution (Density Wave) in the Euler Equations DFR solver and ran convergence studies. The results show a couple of things:
1) Rapid convergence to machine zero across polynomial order and number of elements
2) Order of convergence is N, not the optimal (N+1)
3) Aliasing instability when using the Roe flux
4) The "SlopeMod" limiter destroys the accuracy of the smooth solution

While it was very nice to see (1) happen, the results of (2) made me take a harder look at the DFR approach implemented now. I had taken the "shortcut" of having set the flux at the element faces, and not modifying the interior flux values, which has the impact of converting the flux into a polynomial of order (N-2) due to the removal of the two faces from the order N basis by forcing their values. It seems that the resulting algorithm has a demonstrated order of (N-1) in convergence rate as a result.

One of the reasons I took the above shortcut is that right now we are using the "Legendre-Gauss-Lobato" (LGL) nodes for the algorithm, per Hesthaven, and the LGL points include the face vertices for each element, which makes it impossible to implement the Jameson DFR approach for DFR.

In the Jameson DFR, the Gauss quadrature points do not include the face vertices. The algorithm sets the flux values at the face vertices, then performs an interpolation across the combination of interior and face (N+1+2) vertices to determine the coefficients of the interior flux polynomial such that the new polynomial passes through all (N+1) interior points in the element and the two face vertices. The resulting reconstructed flux polynomial is of order (N+2) and resides on the (N+3) face and solution interior points. Derivatives of this flux are then used directly in the solution algorithm as described in this excellent [AFOSR presentation](http://aero-comlab.stanford.edu/Papers/AFOSR-Meeting-Jul-2014.pdf).

The Jameson DFR algorithm provides an equivalent, but simpler and more efficient (15% faster) way to achieve all of the benefits of DG methods. The DFR solver uses the differential form of the target equations, rather than the integral form, which makes it easier to attack more complex combinations of target equations. DFR has also been extended by Jameson, et al to incorporate spectral methods with entropy stable properties.


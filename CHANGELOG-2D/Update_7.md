## Update: (May 12, 2020): DFR and Aliasing, Instability and Fixing it

On the path to implementing direct flux reconstruction, I found what appeared to be 2nd order aliasing without a clear origin. After consulting Hesthaven(2007) section 5.3, I see that the issue is the interpolation of the flux and subsequently taking the derivative of that interpolated flux. As stated: "the derivative of the interpolation is not the same as the interpolation of the derivative", or put another way, by simply computing the flux from the nodal points of the solution polynomial, we are not treating the flux formally as a polynomial - instead we should perform a formal polynomial fit (projection, instead of interpolation) of the flux prior to using that polynomial to compute derivatives of the flux. The result of using interpolation shown in the text is that we produce an aliasing error into the solution, and their answer is to filter it away instead of using the much more compute intensive projection. The aliasing error also gets worse with increasing N, so the filter should change with N.

I've implemented a simple 2nd order dissipative filter with a constant coefficient, which works well to knock out the oscillations without introducing solution error. I also experimented with a combined 2nd / 4th order filter, similar to what is commonly used in finite volume schemes and for this problem the 2nd order dissipation was enough. However, I can not get stable results beyond about N=6 with this filter, so I'll also investigate an efficient way to do projection and/or improved filters.

--> prior updates
DFR for Maxwell's equations now uses 2nd order artificial dissipation to remove the odd/even aliasing and it works pretty well, which affirms my earlier finding of higher order modes. Next steps might include doing a stability analysis on the DFR scheme, covering the details of the reconstruction. It is slightly different than in Jameson(2014), in that this formulation of NDG uses N+1 points for the polynomials, including the node edges. In Jameson(2014), they extended the flux points to N+3 to cover the element edges, where here I just formed the fluxes at the face points, which are at i=1 and i=N+1. 

You can see it as model 4 with default parameters like ```gocfd -model 4 -graph``` The issue is the instability - it scales with N, so is a higher order modal signal - It could be that the Lax Friedrich's flux, which is 1st order, is not providing damping of the higher order modes. But why are there higher order modes? Are these unresolved waves? I can't answer without more looking into it...



[Back to Index](../CHANGELOG-2D.md)

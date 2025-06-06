## Update: Sep 9 2020):
I've now validated the RT element up to 7th order for divergence of polynomial vector fields. Happily, the special case of zero divergence is being captured with high precision.

The reason this took me quite a while: the choice of basis functions for the RT polynomial must be approached with great care, as the derivatives of the basis must be normalized properly. The literature provides many examples of alternatives for the basis, many of which seem to produce erroneous divergence values that do not show convergence, or even diverge with increasing element order. Others use moments of dot products for the interior, which is not the same process that is used by Jameson and Romero's DFR work. In the interest of staying close to that work, I ended up using the same 2D polynomial basis for the Simplex (triangle) as used by Westhaven with the procedure described by Romero to create the basis for the RT element, which preserves the connection with the Lagrange element interior points. I now have validated that this combination accurately reproduces divergence for polynomial functions and shows convergent behavior with good accuracy for transcendental fields up to 7th order.



[Back to Index](../CHANGELOG-2D.md)

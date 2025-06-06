## Update: Aug 20 2020):
I finally have a working Raviart-Thomas element at any K up to 7th order now. The limitation to 7th order is a matter of
finding optimized point distributions for this element type. This RT element closely follows the approach described in 
[Romero and Jameson](../research/convergence_and_fluxes/DFR/romero_jsc_2017.pdf) for the 3rd and 4th orders, and uses 
the point distributions from
[Williams and Shun](../research/convergence_and_fluxes/DFR/williams-shun-jameson-quadrature.pdf) for orders 1,2 and 5-7.
The RT vector basis cooefficients for order N are solved numerically in one matrix inversion, yielding a full polynomial
basis with (N+1)(N+3) terms implemented as degrees of freedom defining the basis distributed on the [-1,-1] reference 
triangle. There are 3(N+1) points on the faces of the triangle, and N(N+1)/2 points defining the interior, and each of 
the interior points hosts 2 degrees of freedom, one for each orthogonal basis vector [r,s] = [1,0] and [0,1], the axis 
vectos for r and s. The zero order (N=0) element has 0 interior points and 3 face points, one on each face, the N=1 
element has 2 points on each face and 1 point in the interior, the N=2 has 9 points on faces (3 per edge), and 3 interior
points, for a total of (2+1)(2+3)=15 degrees of freedom.

|         RT0 Element         |         RT7 Element         |
|:---------------------------:|:---------------------------:|
| ![](../images/RT1_element.PNG) | ![](../images/RT7_element.PNG) |

The Raviart-Thomas elements at orders 1 and 7 are shown above, with basis vectors as green lines. On the left is RT1, 
showing a single interior point with two basis vectors. The RT7 element has 28 interior points, each with two basis vectors. The edges of RT7 each have 8 points, each with a single normal vector for the basis.

Here is the cool part: the number of interior points within each of the above elements matches that of the Lagrange 
element used for the Nodal Galerkin models(!) That means we can have a Nodal Galerkin element at order N=1, which has 
(N+1)(N+2)/2 = 3 points inside it have a companion RT_2 element that also has 3 interior points. This makes it possible 
to represent gradient, curl, and divergence of our scalar variables residing inside the Nodal Galerkin elements in the 
RT_N+1 elements, where those fields are calculated using the two dimensional polynomial basis in a vector preserving 
way. All this without interpolating between the Lagrangian element and the RT element - the point locations are shared 
between them, which removes a major source of error and complication. 

I'm still implementing tests, but the basics are finally there. Next step: implement gradient, divergence, curl 
operators and test for convergence to known solutions to [ sin(y), sin(x) ]. After that, I'll need to implement the 
calculation of the coordinate transform Jacobian, which actually follows from the previous operators in that the 
calculation of [ J ] includes the same derivatives used in the gradient.



[Back to Index](../CHANGELOG-2D.md)

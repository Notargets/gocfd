# Redefinition of the Raviart Thomas Element Construction (Part 3)

## Issue with the RT3+ element interior basis polynomial & resolution
I've completed an Ervin basis for the RT element to the point where I am 
evaluating the basis matrix and inverse used to project the target flux vector
field onto the basis. There is a problem with RT3+ in that the interior 
basis is not orthogonal, and so the inverted basis projection matrix is 
either singular or ill-conditioned. The issue I've discovered is that the 
2D polynomial field I used to multiply the base interior vectors in the basis 
is the 2D Jacobi polynomial series, parameterized by the Beta parameter. 
That polynomial series is not sufficiently distinct to allow for an 
orthogonal basis.

For RT1, only the two base interior vectors are used. Ervin chose these two 
vectors E4 and E5 for their property that the dot product of each with 
the edge normals on the triangle edges is zero. This ensures that the interior 
basis functions will not directly contribute a flux of the vector field out of 
the element. The Ervin RT1 basis matrix is well conditioned.

For RT2, the E4 and E5 basis vectors are multiplied by polynomials (1-eta-xi),
(xi) and (eta). These effectively produce an orthogonal basis for RT2 and 
the basis matrix is well conditioned.

For the general case of RT3+ (RTk), the E4 and E5 basis vectors are 
multiplied by a 2D polynomial (b)j, which is unique for each interior node 
position. I chose to use 2D Jacobi polynomials, which are parameterized with 
the Beta factor that skews one of the basis polynomial series toward or away 
from the edges. Unfortunately, the RT3+ basis matrices are not well 
conditioned when using these polynomials. After experimentation with various 
Beta parameterization, I was able to get an RT3 basis that was well 
conditioned, but as the element polynomial degree increased for RT4+, the 
conditioning became worse.

The resolution for this is clearly to use a fully orthogonal 2D polynomial - 
the Lagrange polynomial for the simplex. This will extend the same principle 
used in the RT2 basis to RT3+ and lead to well conditioned basis matrices.

## Optimization: Replace edge basis vectors with constant edge normal vectors

The reason Ervin used the three basis vectors E1,E2 and E3 for the edges is 
that they have the property that each of them is independent of the other 
at the midpoint of each edge so that E1 dot E2 is zero at the midpoint of 
edge 2, etc. This property establishes the independence of the basis vectors 
along the edge of the triangle element, specifically for the case of RT0, 
where there are only these basis vectors present.

In the case of RT1+, the edge basis vectors are multiplied by a 
Lagrange polynomial that establishes the same independence between the basis 
vectors, regardless of their form. In other words: for RT1+, we can use any 
edge oriented vectors as the basis, as long as they comply with the requirements
for the RT1+ elements.

The most convenient choice for the edge basis vectors are the edge unit normals.
This allows us to simply project the normal flux component on the basis to 
use the RT1+ elements, avoiding a bunch of extra work to project tangential 
flux onto the edges.

I'll implement a modified version of the RT element that uses unit normal 
vectors and verify proper accuracy in testing.

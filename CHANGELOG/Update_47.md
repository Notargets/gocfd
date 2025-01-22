# Redefinition of the Raviart Thomas Element Construction

I've reviewed the previous construction and it was wrong in substantial ways.
Following is a revised definition of the element, consistent with it's 
purpose and acknowledged attributes.

## Scope of errors in previous declaration(s)

There was a first RT element definition that was very wrong in both how the 
basis was defined, and how the construction was done.

The second iteration of the RT element was correct in the definition of the 
basis, but the construction had a major error in how the evaluation matrix 
was defined. In particular, the basis evaluation matrix was built as a 
Vendermonde matrix, so instead of evaluating the basis function fully at 
each column, each column evaluated the component polynomials of the basis 
function.position

## Revised definition: Jacobi polynomials, evaluated correctly

In this new definition, we'll use the same polynomials used in the scalar 
element, and we'll build it correctly using a full evaluation of the basis 
functions at each position for each DOF.

```
Definition:
	Raviart Thomas element
	Total of (N+1)(N+3) degrees of freedom. (N)(N+1) interior degrees of
	freedom split between R and S directions. 3*(N+1) edge degrees of
	freedom, (N+1) per each of 3 edges on the triangle.

Inputs:
	(N)(N+2)/2 [r,s] points from the interior of the [-1,1] triangle

Outputs:
	[R,S]: Coordinates of points within element

	First (N)(N+1)/2 points: Interior points, excluding edges in
		[-1,1] triangle coordinates corresponding to the R
		direction DOFs
	Next (N)(N+1)/2 points: Same as above for the S direction
	Next (N+1) points: Edge 1 locations
	Next (N+1) points: Edge 2 (Hypotenuse) locations
	Next (N+1) points: Edge 3 locations


Methods:
	The evaluation matrix rt.A relates the coefficients of the RT element
	basis to the solution vectors. For the RT element at degree 1 (RT1),
	this equation is:

⎡ φ₁(P₁)  φ₂(P₁)  φ₃(P₁)  φ₄(P₁)  φ₅(P₁)  φ₆(P₁)  φ₇(P₁)  φ₈(P₁) ⎤   ⎡ c₁ ⎤   ⎡ f₁ ⎤
⎢ φ₁(P₁)  φ₂(P₁)  φ₃(P₁)  φ₄(P₁)  φ₅(P₁)  φ₆(P₁)  φ₇(P₁)  φ₈(P₁) ⎥   ⎢ c₂ ⎥   ⎢ f₂ ⎥
⎢ φ₁(P₂)  φ₂(P₂)  φ₃(P₂)  φ₄(P₂)  φ₅(P₂)  φ₆(P₂)  φ₇(P₂)  φ₈(P₂) ⎥   ⎢ c₃ ⎥   ⎢ f₃ ⎥
⎢ φ₁(P₃)  φ₂(P₃)  φ₃(P₃)  φ₄(P₃)  φ₅(P₃)  φ₆(P₃)  φ₇(P₃)  φ₈(P₃) ⎥   ⎢ c₄ ⎥   ⎢ f₄ ⎥
⎢ φ₁(P₄)  φ₂(P₄)  φ₃(P₄)  φ₄(P₄)  φ₅(P₄)  φ₆(P₄)  φ₇(P₄)  φ₈(P₄) ⎥ * ⎢ c₅ ⎥ = ⎢ f₅ ⎥
⎢ φ₁(P₅)  φ₂(P₅)  φ₃(P₅)  φ₄(P₅)  φ₅(P₅)  φ₆(P₅)  φ₇(P₅)  φ₈(P₅) ⎥   ⎢ c₆ ⎥   ⎢ f₆ ⎥
⎢ φ₁(P₆)  φ₂(P₆)  φ₃(P₆)  φ₄(P₆)  φ₅(P₆)  φ₆(P₆)  φ₇(P₆)  φ₈(P₆) ⎥   ⎢ c₇ ⎥   ⎢ f₇ ⎥
⎣ φ₁(P₇)  φ₂(P₇)  φ₃(P₇)  φ₄(P₇)  φ₅(P₇)  φ₆(P₇)  φ₇(P₇)  φ₈(P₇) ⎦   ⎢ c₈ ⎥   ⎢ f₈ ⎥

Notes for RT1:
	1) there is one internal position in [R,S] shared by the two interior DOFs,
	one for each of R and S directions.
	2) there are 6 positions, one for each of the DOFs on the edges, 2
	per edge.
	3) each basis function φ is evaluated fully at each position. This
	amounts to summation of the full polynomial range of each function
	at the corresponding position.

	Note that (3) above is in contrast to the Vendermonde matrix, which
	is very different in that each column in a Vandermonde matrix is an
	evaluation of a single polynomial term. The Vandermonde matrix
	relates the coefficients of the polynomial terms to the function
	representation in a polynomial, which is NOT what we are doing here.

	For an RT element, the function is distinct on the interior and the
	edges. There is no value of the edge basis functions in the
	interior, and there is no value of the interior functions on the
	edge. We correct the above evaluation matrix for RT1 to reflect the
	nature of the RT element:

⎡ φ₁(P₁)  φ₂(P₁)    0       0       0       0       0       0    ⎤   ⎡ c₁ ⎤   ⎡ f₁ ⎤
⎢ φ₁(P₁)  φ₂(P₁)    0       0       0       0       0       0    ⎥   ⎢ c₂ ⎥   ⎢ f₂ ⎥
⎢   0       0     φ₃(P₂)  φ₄(P₂)    0       0       0       0    ⎥   ⎢ c₃ ⎥   ⎢ f₃ ⎥
⎢   0       0     φ₃(P₃)  φ₄(P₃)    0       0       0       0    ⎥   ⎢ c₄ ⎥   ⎢ f₄ ⎥
⎢   0       0       0       0     φ₅(P₄)  φ₆(P₄)    0       0    ⎥   ⎢ c₅ ⎥ = ⎢ f₅ ⎥
⎢   0       0       0       0     φ₅(P₅)  φ₆(P₅)    0       0    ⎥   ⎢ c₆ ⎥   ⎢ f₆ ⎥
⎢   0       0       0       0       0       0     φ₇(P₆)  φ₈(P₆) ⎥   ⎢ c₇ ⎥   ⎢ f₇ ⎥
⎢   0       0       0       0       0       0     φ₇(P₇)  φ₈(P₇) ⎥   ⎢ c₈ ⎥   ⎢ f₈ ⎥

Here is the RT2 evaluation matrix to make the pattern clear.
Notes for RT2:
	1) There are 3 internal positions in [R,S] (N=2, (N)*(N+1)/2 = 3
	2) There are 3 edge DOFs for each edge, total of 9 positions
	3) The edge basis functions are only evaluated on each edge and the interior
	basis functions are only evaluated at each interior position.

⎡ φ₁(P₁)  φ₂(P₁)  φ₃(P₁)  φ₄(P₁)  φ₅(P₁)  φ₆(P₁)    0       0       0       0       0       0    0        0        0     ⎤   ⎡ c₁  ⎤   ⎡ f₁  ⎤
⎢ φ₁(P₂)  φ₂(P₂)  φ₃(P₂)  φ₄(P₂)  φ₅(P₂)  φ₆(P₂)    0       0       0       0       0       0    0        0        0     ⎥   ⎢ c₂  ⎥   ⎢ f₂  ⎥
⎢ φ₁(P₃)  φ₂(P₃)  φ₃(P₃)  φ₄(P₃)  φ₅(P₃)  φ₆(P₃)    0       0       0       0       0       0    0        0        0     ⎥   ⎢ c₃  ⎥   ⎢ f₃  ⎥
⎢ φ₁(P₁)  φ₂(P₁)  φ₃(P₁)  φ₄(P₁)  φ₅(P₁)  φ₆(P₁)    0       0       0       0       0       0    0        0        0     ⎥   ⎢ c₄  ⎥   ⎢ f₄  ⎥
⎢ φ₁(P₂)  φ₂(P₂)  φ₃(P₂)  φ₄(P₂)  φ₅(P₂)  φ₆(P₂)    0       0       0       0       0       0    0        0        0     ⎥   ⎢ c₅  ⎥   ⎢ f₅  ⎥
⎢ φ₁(P₃)  φ₂(P₃)  φ₃(P₃)  φ₄(P₃)  φ₅(P₃)  φ₆(P₃)    0       0       0       0       0       0    0        0        0     ⎥   ⎢ c₆  ⎥   ⎢ f₆  ⎥
⎢   0       0       0       0       0       0     φ₇(P₄)  φ₈(P₄)  φ₉(P₄)    0       0       0    0        0        0     ⎥   ⎢ c₇  ⎥   ⎢ f₇  ⎥
⎢   0       0       0       0       0       0     φ₇(P₅)  φ₈(P₅)  φ₉(P₅)    0       0       0    0        0        0     ⎥   ⎢ c₈  ⎥ = ⎢ f₈  ⎥
⎢   0       0       0       0       0       0     φ₇(P₆)  φ₈(P₆)  φ₉(P₆)    0       0       0    0        0        0     ⎥   ⎢ c₉  ⎥   ⎢ f₉  ⎥
⎢   0       0       0       0       0       0       0       0       0 φ₁₀(P₇) φ₁₁(P₇) φ₁₂(P₇)    0        0        0     ⎥   ⎢ c₁₀ ⎥   ⎢ f₁₀ ⎥
⎢   0       0       0       0       0       0       0       0       0 φ₁₀(P₈) φ₁₁(P₈) φ₁₂(P₈)    0        0        0     ⎥   ⎢ c₁₁ ⎥   ⎢ f₁₁ ⎥
⎢   0       0       0       0       0       0       0       0       0 φ₁₀(P₉) φ₁₁(P₉) φ₁₂(P₉)    0        0        0     ⎥   ⎢ c₁₂ ⎥   ⎢ f₁₂ ⎥
⎢   0       0       0       0       0       0       0       0       0       0       0       0  φ₁₃(P₁₀) φ₁₄(P₁₀) φ₁₅(P₁₀)⎥   ⎢ c₁₃ ⎥   ⎢ f₁₃ ⎥
⎢   0       0       0       0       0       0       0       0       0       0       0       0  φ₁₃(P₁₁) φ₁₄(P₁₁) φ₁₅(P₁₁)⎥   ⎢ c₁₄ ⎥   ⎢ f₁₄ ⎥
⎢   0       0       0       0       0       0       0       0       0       0       0       0  φ₁₃(P₁₂) φ₁₄(P₁₂) φ₁₅(P₁₂)⎥   ⎢ c₁₅ ⎥   ⎢ f₁₅ ⎥
```
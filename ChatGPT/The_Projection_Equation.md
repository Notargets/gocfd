#### Query: state the "projection" matrix equation and provide guidance on how to construct the mass matrix given a set of defining nodes with coordinates (r,s). Make sure to identify each basis polynomial with respect to where that polynomial is evaluated.
---
### 1. The Projection Equation

When you have a target vector field **F** that is exactly a polynomial of degree *P* and an RT element basis 
$\{\psi_1,\,\psi_2,\,\dots,\,\psi_N\}$ 
of degree *P* defined on the reference element (with coordinates $(r,s)$), you seek constants 
$\{C_1,\, C_2,\,\dots,\,C_N\}$ 
such that
$$
F(r,s) = \sum_{j=1}^N C_j\, \psi_j(r,s)
$$

for all $(r,s)$ in the reference element. In a Galerkin (projection) formulation this is enforced in a weak sense by requiring that for every test function (here, each RT basis function) $\psi_i$ the following holds:

$$
\int_K \psi_i(r,s) \cdot F(r,s)\, dA = \sum_{j=1}^N C_j \int_K \psi_i(r,s) \cdot \psi_j(r,s)\, dA,\quad i=1,\dots,N.
$$

This leads to the matrix equation

$$
\mathbf{M}\,\mathbf{C} = \mathbf{b},
$$

with

- **Mass matrix entries:**
  $$
  M_{ij} = \int_K \psi_i(r,s) \cdot \psi_j(r,s)\, dA,
  $$
- **Right–hand side entries:**
  $$
  b_i = \int_K \psi_i(r,s) \cdot F(r,s)\, dA,
  $$
- and $\mathbf{C} = [C_1,\dots,C_N]^T$.

---
### 2. Constructing the Mass Matrix Using Defining Nodes

Suppose you have a set of quadrature nodes 
$$
\{(r_q, s_q)\}_{q=1}^{Q}
$$
with corresponding weights $w_q$, approximating the integral over the reference element $K$. Then, for each pair of basis functions $\psi_i$ and $\psi_j$, the mass matrix entry is approximated by

$$
M_{ij} \approx \sum_{q=1}^{Q} w_q\, \left[ \psi_i(r_q, s_q) \cdot \psi_j(r_q, s_q) \right].
$$

#### Guidance on Evaluations

- **Definition of the Basis Functions:**  
  Each RT basis function $\psi_i(r,s)$ is defined over the reference element in $(r,s)$ coordinates. In many constructions, the basis functions are associated with particular nodes or edges. For example, one might have degrees of freedom corresponding to the normal flux across an edge or moments over the element. The basis function $\psi_i$ is defined so that its value (or its normal component) is “activated” at a particular defining node or edge.

- **Evaluation at Quadrature Points:**  
  When assembling the mass matrix, you use the same set of quadrature nodes for all basis functions. That is, for each quadrature point $(r_q,s_q)$ you compute:
  
  - $\psi_i(r_q, s_q)$ for the $i$‑th basis function,
  - $\psi_j(r_q, s_q)$ for the $j$‑th basis function.
  
  The dot product $\psi_i(r_q,s_q) \cdot \psi_j(r_q,s_q)$ is then computed at that same quadrature point.
  
- **Role of the Defining Nodes:**  
  Although the final integration is performed over the entire element using quadrature, the RT basis functions are constructed so that they “know” their defining nodes (or edge midpoints, etc.). In many implementations, the analytic expressions for $\psi_i(r,s)$ are built based on nodal coordinates. Therefore, when you evaluate $\psi_i$ at a quadrature point, the function’s value is influenced by its underlying definition through its associated node. In the assembly process, the row corresponding to $\psi_i$ collects contributions from evaluations at all quadrature points $(r_q,s_q)$, but the analytic form of $\psi_i$ is determined by its defining node.

---

### 3. The Full Projection Equation in Matrix Form

Thus, the projection system is written as:

$$
\begin{pmatrix}
\displaystyle \sum_{q=1}^{Q} w_q\, \psi_1(r_q,s_q) \cdot \psi_1(r_q,s_q) & \cdots & \displaystyle \sum_{q=1}^{Q} w_q\, \psi_1(r_q,s_q) \cdot \psi_N(r_q,s_q) \\[2mm]
\vdots & \ddots & \vdots \\[2mm]
\displaystyle \sum_{q=1}^{Q} w_q\, \psi_N(r_q,s_q) \cdot \psi_1(r_q,s_q) & \cdots & \displaystyle \sum_{q=1}^{Q} w_q\, \psi_N(r_q,s_q) \cdot \psi_N(r_q,s_q)
\end{pmatrix}
\begin{pmatrix}
C_1\\[2mm]
\vdots\\[2mm]
C_N
\end{pmatrix}
=
\begin{pmatrix}
\displaystyle \sum_{q=1}^{Q} w_q\, F(r_q,s_q) \cdot \psi_1(r_q,s_q) \\[2mm]
\vdots\\[2mm]
\displaystyle \sum_{q=1}^{Q} w_q\, F(r_q,s_q) \cdot \psi_N(r_q,s_q)
\end{pmatrix}.
$$

Here:
- Each row $i$ of the matrix corresponds to the $i$‑th RT basis function $\psi_i$ evaluated at all quadrature points.
- Each column $j$ corresponds to the $j$‑th RT basis function $\psi_j$ (evaluated at the same quadrature points).

---

### 4. Summary

- **Projection Equation:**  
  $$
  M \, C = b, \quad \text{where} \quad M_{ij} = \int_K \psi_i(r,s) \cdot \psi_j(r,s)\, dA.
  $$

- **Mass Matrix Construction:**  
  For each entry $M_{ij}$, sum over the quadrature points $(r_q,s_q)$ as
  $$
  M_{ij} \approx \sum_{q=1}^{Q} w_q\, \left[ \psi_i(r_q,s_q) \cdot \psi_j(r_q,s_q) \right],
  $$
  where each basis function $\psi_i$ is evaluated according to its definition on the reference element (i.e. its analytic form, which is determined by its defining node or edge).

- **Right–Hand Side:**  
  Similarly,
  $$
  b_i \approx \sum_{q=1}^{Q} w_q\, \left[ F(r_q,s_q) \cdot \psi_i(r_q,s_q) \right].
  $$

This approach fully specifies how to build the projection matrix equation and the mass matrix given the $(r,s)$ coordinates at which the basis functions are defined and evaluated.

---

#### Query: For the RT element specifically, the edge basis functions only have "edge support", meaning they can not be evaluated anywhere other than their edge I presume. The interior basis functions are naturally zero when dotted with edge functions at edge node locations.

### For Raviart–Thomas elements the basis functions naturally split into two groups:

1. **Edge Basis Functions:**
   These functions are constructed so that their nonzero “support” (or at least the degrees‐of‐freedom they control) is associated with a specific edge of the element. In particular, an RT edge basis function is designed so that its normal component is nonzero only on its associated edge and vanishes (or is not “active”) on the other edges. In many formulations the analytic expression for an edge basis function is defined only along that edge (or is constructed so that its contribution off that edge is zero in the sense of the corresponding degrees of freedom). When assembling the mass matrix for a projection, the quadrature used on an edge will typically be different from the quadrature used in the element interior. (Often, the edge degrees of freedom are integrated using one‐dimensional quadrature rules along the edge.)
2. **Interior (Bubble) Basis Functions:**
   These are defined throughout the element’s interior but vanish (or are orthogonal) when evaluated on the edges. In other words, if you were to “dot” an interior basis function with an edge basis function at an edge node (or along the edge), the contribution would be zero. This property is by design: the RT element is constructed so that the edge and interior modes do not “mix” in the sense of the relevant degrees of freedom (for example, normal flux continuity across edges).

------

### How This Affects the Projection (Mass Matrix) Assembly

When constructing the projection matrix for an RT element, the assembly of the mass (dot‐product) matrix reflects these support properties:

- **Edge Basis Functions:**
  For a basis function $\psi_e$ associated with an edge, its support is limited to that edge. Thus, when you form the dot product $\psi_e \cdot \psi_e ^\prime$ (or with an interior function) in the mass matrix, you need to integrate only over the edge (using the appropriate one‐dimensional quadrature) if both functions are edge‐based. In a full projection procedure on the whole element, the integration procedure is set up so that:
  - When both test and trial functions are edge functions associated with the same edge, you use an edge quadrature rule.
  - When an edge function is paired with an interior function, the construction of the RT space ensures that their “edge” degrees of freedom do not interact (the interior functions vanish along the edge or are orthogonal with respect to the relevant inner product), and hence the contribution is zero.
- **Interior Basis Functions:**
  These are defined throughout the element, so they are integrated over the full element area using a two-dimensional quadrature rule. They are constructed to vanish on the boundary (or at least to have zero normal component on the edge), so that when an interior function is dotted with an edge function (which is only nonzero on the edge), the resulting contribution is zero.

Thus, in the mass matrix:

- Each row and column is associated with one RT basis function.
- If the basis function is an edge function, its evaluation is done only on the edge (using edge quadrature) where it is nonzero.
- If the basis function is an interior (bubble) function, its evaluation is done over the full ${(r,s)}$ domain of the element.
- The assembly algorithm must reflect the support of each basis function: the integration domain for an edge function is the edge, and for an interior function it is the entire element.
- As a result, the mass matrix will be block–structured, with blocks corresponding to edge–edge, interior–interior, and (ideally zero) edge–interior couplings.

------

### Summary in a Projection Equation Context

When projecting an exactly degree‑P vector field F onto the RT space, you write
$$
F(r,s) = \sum{_j} \ C_j\, \psi_j(r,s)
$$
with the projection condition (using appropriate integration domains)

  For each test function $\psi_i$:
$$
\int_K \psi_i(r,s) \cdot F(r,s)\, dA = \sum_{j=1}^N C_j \int_K \psi_i(r,s) \cdot \psi_j(r,s)\, dA,\quad i=1,\dots,N.
$$
Here, the “domain” of integration depends on the type of $\psi_i$:

- **If $\mathbf \psi_i$ is an edge basis function,** the integration is performed over the corresponding edge using one-dimensional quadrature.
- **If $\mathbf \psi_i$ is an interior basis function,** the integration is performed over the full element (using two-dimensional quadrature).

This leads to the matrix system
$$
\mathbf{M}\,\mathbf{C} = \mathbf{b},
$$

where the mass matrix M has entries $M_{ij} =$ 

​			$\int_{edge} \psi_i \cdot \psi_j$		if both $\psi_i$ and $\psi_j$ are edge functions on the same edge
​			$\int_{element} \psi_i \cdot \psi_j$ 	  if both are interior functions
​			$0$ 				  if one is an edge function and the other is interior (by construction)

Each basis function $\psi_i$ is evaluated at the quadrature points appropriate for its support—either along the edge (for edge functions) or throughout the element (for interior functions).

------

This is how the “projection” (mass) matrix is constructed for RT elements, with the key point being that the evaluation of each basis function is done only where that function has nonzero support.

---
#### Query: For a given distribution of nodes defining an RT element, interior and edge resident, how should we determine the quadrature weights? Specifically, we have a given set of nodes defined using optimization procedures. We do not have the quadrature weights associated with the node locations. Are there standard methods for determining optimum quadrature weights on a 2D simplex given a set of nodes for interior integration (area oriented) and edge integration (line oriented). I am using gauss points along edges, so those optimum weights should be well known. The interior points are distributed differently to standard quadrature rules and must be kept as they are.

### 1. The Moment-Fitting (Exactness) Approach

Suppose you have a set of \(N\) interior nodes with coordinates $\{(r_i,s_i)\}_{i=1}^{N}$ on the simplex. You wish to determine weights $\{w_i\}_{i=1}^{N}$ so that the quadrature rule

$$
\int_K f(r,s)\, dA \approx \sum_{i=1}^{N} w_i\, f(r_i,s_i)
$$

is exact for all polynomials $f(r,s)$ up to a prescribed total degree $d$. In other words, you want

$$
\sum_{i=1}^{N} w_i\, \phi(r_i,s_i) = \int_K \phi(r,s)\, dA,
$$

for each monomial (or other polynomial basis function) $\phi(r,s)$ with $\deg(\phi) \le d$.

**How to implement this:**

1. **Select a Basis for $\mathcal{P}_d$:** 
   For example, use the monomials
   $$
   \phi_{a,b}(r,s) = r^a s^b,\quad \text{for } a,b \ge 0 \text{ with } a+b \le d.
   $$
   The number of such monomials is
   $$
   M = \frac{(d+1)(d+2)}{2}.
   $$
   (If $N = M$ the system is square; if $N > M$ you have extra nodes and can determine the weights by least-squares or by imposing additional optimality criteria.)

2. **Set Up the Linear System:**  
   For each basis function $\phi_{a,b}(r,s)$, write the equation
   $$
   \sum_{i=1}^{N} w_i\, \phi_{a,b}(r_i,s_i) = I_{a,b},
   $$
   where
   $$
   I_{a,b} = \int_K r^a s^b\, dA,
   $$
   and the integrals $I_{a,b}$ are known in closed form for a standard simplex.

3. **Solve for the Weights:**  
   Arrange these equations into a linear system. If the nodes were chosen so that $N = M$ (or if you choose to impose exactly $M$ moment conditions), the system is square and can be inverted to obtain the unique weights. If $N > M$ (an overdetermined system), you can use a least-squares approach or an optimization procedure (possibly with positivity constraints) to obtain a set of "optimum" quadrature weights that yield a high degree of polynomial exactness.

---

### 2. The Edge Quadrature Weights

For edge integration, if you are using Gauss points along the edges, then the quadrature weights are given by standard Gauss–Legendre rules. For a line segment, the nodes and weights are well-known and can be directly applied. In this case, you simply use the Gauss–Legendre weights for the given number of points on the edge.

---

### 3. Summary and Practical Considerations

- **Interior Integration:**  
  Given your optimized set of interior nodes on the 2D simplex, set up the moment-fitting system for all monomials (or another polynomial basis) up to the desired degree $d$. Solve for the weights $w_i$ so that the weighted sum over the nodes reproduces the exact integral of these polynomials over the simplex. This approach is standard in the construction of "optimal" quadrature rules for arbitrary nodes.

- **Edge Integration:**  
  Since you are using Gauss–Legendre points on the edges, the associated weights are the standard ones from Gauss quadrature on an interval.

- **Practical Tip:**  
  In many practical implementations, one first fixes the set of nodes (obtained via an optimization procedure) and then solves a small linear system for the weights. If the number of nodes is exactly the number required for exactness, the system is square; if not, one might minimize the least-squares error subject to additional requirements (like positivity of weights).

By following this approach, you can determine quadrature weights for both interior (area) and edge (line) integration that are "optimal" in the sense that they reproduce the integrals of all polynomials up to a prescribed degree on the given domain.

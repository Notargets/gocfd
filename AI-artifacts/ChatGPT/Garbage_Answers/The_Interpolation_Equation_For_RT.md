# Computing the Divergence of the Vector Field $\mathbf{F}$

## 1. Structure of the Basis Functions
We express the vector field $\mathbf{F}$ in terms of the Raviart-Thomas basis functions:

$$
\mathbf{F}(r,s) = \sum_{j=1}^{N_{\text{int}}} c_j^{(E4)} \boldsymbol{\psi}_j^{(E4)}(r,s) + c_j^{(E5)} \boldsymbol{\psi}_j^{(E5)}(r,s)
+ \sum_{k=1}^{N_{\text{edge}}} c_k^{(n)} \boldsymbol{\psi}_k^{(n)}(r,s)
$$

where:
- $N_{\text{int}} = \frac{(P+1)(P+2)}{2}$ is the number of interior DOFs (two per interior node).
- $N_{\text{edge}} = 3(P+1)$ is the number of edge DOFs (one per edge node for normal continuity).

The total number of DOFs is:

$$
N = N_{\text{int}} + N_{\text{edge}} = (P+1)(P+3)
$$

## 2. Interpolation Equations
Each DOF provides a linear constraint, forming a system of equations.

### (a) Interior Nodes
For each interior node $(r_i, s_i)$, we impose:

$$
\mathbf{F}(r_i, s_i) \cdot \mathbf{E4}(r_i, s_i) = d_i^{(E4)}
$$

$$
\mathbf{F}(r_i, s_i) \cdot \mathbf{E5}(r_i, s_i) = d_i^{(E5)}
$$

Expanding using only interior basis functions:

$$
\sum_{j=1}^{N_{\text{int}}} c_j^{(E4)} \boldsymbol{\psi}_j^{(E4)}(r_i, s_i) \cdot \mathbf{E4}(r_i, s_i)
+ \sum_{j=1}^{N_{\text{int}}} c_j^{(E5)} \boldsymbol{\psi}_j^{(E5)}(r_i, s_i) \cdot \mathbf{E5}(r_i, s_i) = d_i^{(E4)}
$$

$$
\sum_{j=1}^{N_{\text{int}}} c_j^{(E4)} \boldsymbol{\psi}_j^{(E4)}(r_i, s_i) \cdot \mathbf{E5}(r_i, s_i)
+ \sum_{j=1}^{N_{\text{int}}} c_j^{(E5)} \boldsymbol{\psi}_j^{(E5)}(r_i, s_i) \cdot \mathbf{E5}(r_i, s_i) = d_i^{(E5)}
$$

### (b) Edge Nodes
At each edge node $(r_k, s_k)$, the normal component is enforced:

$$
\mathbf{F}(r_k, s_k) \cdot \mathbf{n}_k = g_k
$$

Expanding using only edge basis functions:

$$
\sum_{j=1}^{N_{\text{edge}}} c_j^{(n_k)} \boldsymbol{\psi}_j^{(n_k)}(r_k, s_k) \cdot \mathbf{n}_k = g_k
$$

### (c) Zero Cross-Terms
Since interior basis functions $\boldsymbol{\psi}_j^{(E4)}$, $\boldsymbol{\psi}_j^{(E5)}$ are constructed only within the element interior, they satisfy:

$$
\boldsymbol{\psi}_j^{(E4)} \cdot \mathbf{n}_k = 0, \quad
\boldsymbol{\psi}_j^{(E5)} \cdot \mathbf{n}_k = 0
$$

This means that the edge constraints only involve edge basis functions.

## 3. Matrix System
Now, we can express the full system in block matrix form:

$$
\begin{bmatrix}
A_{\text{int}} & 0 \\
0 & A_{\text{edge}}
\end{bmatrix}
\begin{bmatrix}
\mathbf{c}_{\text{int}} \\
\mathbf{c}_{\text{edge}}
\end{bmatrix}
=
\begin{bmatrix}
\mathbf{b}_{\text{int}} \\
\mathbf{b}_{\text{edge}}
\end{bmatrix}
$$

where:
- $A_{\text{int}}$ represents interior projection equations (purely involving $\mathbf{E4}$, $\mathbf{E5}$).
- $A_{\text{edge}}$ represents edge normal constraints (purely involving $\mathbf{n}$).
- The block structure arises because interior and edge basis functions are independent.

The unknown coefficients are split:

$$
\mathbf{c} =
\begin{bmatrix}
\mathbf{c}_{\text{int}} \\
\mathbf{c}_{\text{edge}}
\end{bmatrix}
$$

The right-hand side is:

$$
\mathbf{b} =
\begin{bmatrix}
\mathbf{b}_{\text{int}} \\
\mathbf{b}_{\text{edge}}
\end{bmatrix}
$$

## 4. Computing the Divergence
To compute the divergence of $\mathbf{F}$:

$$
\nabla \cdot \mathbf{F}(r,s) = \sum_{j=1}^{N} c_j \nabla \cdot \boldsymbol{\psi}_j(r,s)
$$

This can be represented in matrix form:

$$
\text{div}(\mathbf{F}) = D \mathbf{c}
$$

where:
- $D$ is the divergence matrix with entries $D_{ij} = \nabla \cdot \boldsymbol{\psi}_j(r_i, s_i)$.
- $\mathbf{c}$ is the vector of coefficients combining interior and edge contributions.

Since the divergence matrix $D$ only depends on the basis functions, it can be precomputed, and the divergence field at specific points is obtained through this matrix-vector multiplication.

I’m considering how to construct the test shock wave - perhaps I can use the analytical shock definition on the RT element with the edge nodes, etc, then attenuate the modal coefficients of a Vandermonde matrix for the RT element and re-project the nodal values using the RT matrix for interpolation to all RT nodes?

---



# Steps for Constructing an RT-Compatible Shock Profile

## Step 1: Define the Analytical Shock Profile

We start with a standard smooth shock transition, modeled as a hyperbolic tangent function:

$$
\rho(x) = \rho_1 + \frac{\rho_2 - \rho_1}{2} \left(1 + \tanh\left(\frac{x - x_s}{\epsilon}\right)\right)
$$

where:
- $$ x_s $$ is the shock position,
- $$ \epsilon $$ controls the smoothness of the transition.

For a shock positioned diagonally across the RT element, define:

$$
\rho(r,s) = \rho_1 + \frac{\rho_2 - \rho_1}{2} \left(1 + \tanh\left(\frac{r+s}{\epsilon}\right)\right)
$$

which places the shock along the diagonal $$ r + s = 0 $$ of the reference triangle.

## Step 2: Evaluate at RT Nodes

Using the RT element node set, evaluate the shock profile at each node:

$$
\mathbf{u}_{\text{nodes}} = \left[ u(\mathbf{x}_1), u(\mathbf{x}_2), \dots, u(\mathbf{x}_N) \right]
$$

where $$ u(r,s) $$ represents the shock profile evaluated at each RT node.

## Step 3: Compute Modal Coefficients (Vandermonde Projection)

Project the nodal values into the RT modal basis using the inverse Vandermonde matrix:

$$
\mathbf{c} = V^{-1} \mathbf{u}_{\text{nodes}}
$$

where:
- $$ V $$ is the Vandermonde matrix for the RT basis.
- $$ \mathbf{c} $$ contains the modal coefficients.

## Step 4: Apply Modal Filtering (Persson Dissipation)

Apply an artificial dissipation filter to damp high-order modal coefficients:

$$
c_k \leftarrow \left(1 - \nu \left(\frac{k}{N}\right)^{2p} \right) c_k
$$

where:
- $$ k $$ indexes the modal "degree," not the node number.
- $$ N $$ is the max polynomial degree in the basis.
- $$ \nu $$ is a user-defined dissipation strength (tunable),
- $$ p $$ controls sharpness (typically 1 to 4 — sharper if you only want to hit the highest modes).

For RT, you may want to split this by **type**:

- Edge modes may need lighter filtering (they interact with neighbors directly).
- Interior modes could get heavier filtering.

This is flexible, but the basic form holds.

## Step 5: Reconstruct Nodal Values

After filtering, reconstruct the nodal values by applying the Vandermonde matrix:

$$
\mathbf{u}_{\text{nodes}}^{\text{filtered}} = V \mathbf{c}_{\text{filtered}}
$$

## Step 6: (Optional) Interpolate to Quadrature Points

If you want to check how this looks at quadrature points (for flux projection, edge flux reconstruction, etc.), you can do:

$$
\mathbf{u}_{\text{quad}} = V_{\text{quad}} \mathbf{c}_{\text{filtered}}
$$

where $$ V_{\text{quad}} $$ is the Vandermonde matrix evaluated at quadrature points.

## Summary Table

| Step | Description                                            |
| ---- | ------------------------------------------------------ |
| 1    | Define shock profile in $$ (r,s) $$                    |
| 2    | Evaluate at RT nodes                                   |
| 3    | Project to RT modal coefficients (Vandermonde inverse) |
| 4    | Apply Persson modal filter to coefficients             |
| 5    | Reconstruct filtered nodal values                      |
| 6    | (Optional) Interpolate to quadrature points            |
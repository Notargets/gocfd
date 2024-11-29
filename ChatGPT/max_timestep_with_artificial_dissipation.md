---
title: Discussion with ChatGPT 4o about the maximum stable timestep under a CFL constraint for Persson's artificial viscocity
---

# Question

What is the max eigenvalue of persson’s artificial viscosity for simplex discontinuous galerkin finite elements used to solve the euler equations?

# Answer

The maximum eigenvalue of Persson's artificial viscosity for simplex discontinuous Galerkin finite elements used to solve the Euler equations depends on the specific formulation and normalization of the viscosity model. However, in the commonly used framework proposed by Persson and Peraire (2006), the artificial viscosity $\nu$ is scaled to stabilize high-frequency oscillations while maintaining consistency and conservation properties.

### Key Points:
1. **Definition and Scaling**:
   Persson's artificial viscosity is usually defined locally based on an indicator function $\eta$, which measures the smoothness of the solution. A common form is:
   $$
   \nu_K = C \, h_K \, \lambda_{\text{max}}(u)
   $$
   where:
   - $C$ is a tunable constant,
   - $h_K$ is a characteristic length scale of the element $K$,
   - $\lambda_{\text{max}}(u)$ is the maximum eigenvalue of the physical flux Jacobian.

2. **Flux Jacobian and Maximum Eigenvalue**:
   For the Euler equations in $d$-dimensions, the maximum eigenvalue $\lambda_{\text{max}}$ of the flux Jacobian corresponds to the speed of sound plus the maximum velocity component, which can be expressed as:
   $$
   \lambda_{\text{max}}(u) = |\vec{u}| + c,
   $$
   where:
   - $|\vec{u}|$ is the magnitude of the flow velocity,
   - $c = \sqrt{\gamma p / \rho}$ is the speed of sound,
   - $\gamma$ is the specific heat ratio,
   - $p$ and $\rho$ are the pressure and density.

3. **Maximum Artificial Viscosity**:
   In the simplex DG context:
   $$
   \nu_{\text{max}} = C \, h_{\text{max}} \, (\max |\vec{u}| + \max c),
   $$
   where $h_{\text{max}}$ is the largest element size in the mesh.

The value of $C$ is problem-dependent and typically chosen empirically (e.g., between $0.1$ and $1.0$) to ensure stabilization without overly diffusing sharp features like shocks or contact discontinuities.

If you'd like, I can help derive specific values or expressions for a given problem setup or implementation details!

# Question

I need the maximum eigenvalue for use in deriving the maximum stable timestep within a CFL limit for the time stepping method.

# Answer

To derive the maximum stable timestep within a CFL (Courant-Friedrichs-Lewy) limit for the Euler equations using Persson's artificial viscosity and simplex discontinuous Galerkin finite elements, we must consider the **maximum eigenvalue** of the combined system (physical fluxes + artificial viscosity). Here's the breakdown:

---

### **Physical Contribution: Maximum Eigenvalue**
For the Euler equations in $d$-dimensions:
$$
\lambda_{\text{physical, max}} = |\vec{u}| + c,
$$
where:
- $|\vec{u}|$ is the velocity magnitude ($|\vec{u}| = \sqrt{u_x^2 + u_y^2 + u_z^2}$ in 3D),
- $c = \sqrt{\gamma p / \rho}$ is the speed of sound,
- $\gamma$ is the ratio of specific heats.

---

### **Artificial Viscosity Contribution**
In Persson's artificial viscosity framework, the artificial diffusion term adds an eigenvalue contribution proportional to the element size and the maximum eigenvalue of the physical flux:
$$
\lambda_{\text{viscosity, max}} = \frac{\nu}{h_K},
$$
where $\nu$ is the artificial viscosity coefficient:
$$
\nu = C \, h_K \, \lambda_{\text{physical, max}}.
$$
Substituting this into $\lambda_{\text{viscosity, max}}$, we get:
$$
\lambda_{\text{viscosity, max}} = C \, \lambda_{\text{physical, max}}.
$$

---

### **Combined Maximum Eigenvalue**
The total maximum eigenvalue for the system is the sum of the physical and artificial viscosity contributions:
$$
\lambda_{\text{max, total}} = \lambda_{\text{physical, max}} + \lambda_{\text{viscosity, max}}.
$$
Substitute the expressions for each term:
$$
\lambda_{\text{max, total}} = (1 + C) \, \lambda_{\text{physical, max}}.
$$

Thus, the total maximum eigenvalue is:
$$
\lambda_{\text{max, total}} = (1 + C)(|\vec{u}| + c).
$$

---

### **Maximum Stable Timestep**
The CFL condition imposes the following restriction on the timestep $\Delta t$:
$$
\Delta t \leq \frac{\text{CFL} \, h_K}{\lambda_{\text{max, total}}},
$$
where $\text{CFL}$ is the chosen Courant number (typically $\text{CFL} \leq 1$ for stability).

Substituting $\lambda_{\text{max, total}}$:
$$
\Delta t \leq \frac{\text{CFL} \, h_K}{(1 + C)(|\vec{u}| + c)}.
$$

---

### Final Expression for Maximum Stable Timestep:
$$
\Delta t \leq \frac{\text{CFL} \, h_K}{(1 + C)(|\vec{u}| + \sqrt{\gamma p / \rho})}.
$$

This expression provides the maximum timestep based on the combined contributions of the physical fluxes and Persson's artificial viscosity.

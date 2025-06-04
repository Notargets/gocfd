Isn't it true that if I subdivide the global element edge into P+1 pieces, do a first order flux reconstruction on each of those pieces to get each sub-face value, then re-assemble the pieces to compose each edge node flux value that I have a P-th order reconstruction?

The only information needed for this approach is what I have: the interior node on either side of the piece of the global edge for each sub-face - I have that from the sub-element triangulation. The global edges are formed using Gauss points on the linear edges, so I have the weights needed to do a high order integration. I don't need the node values on the global element vertices or on the shared faces, because that's not needed for the first order flux (or WENO) reconstruction in FV methods. I would like to use the existing interior node locations, which are not equispaced, etc. For a first order or WENO flux reconstruction, how would I normally proceed with a FV based approach, given the following: - The single interior nodes opposing the sub-element face; one from the neighbor element interior and one from the current element interior - For WENO, other information from interior nodes to each global element, as needed.

# **1. First-Order (Godunov-Type) Flux Reconstruction**
### **Step 1: Identify the Two Opposing Interior Nodes for Each Sub-Face**
For each sub-face of the global edge, determine:
- **$$U_L$$** → Value from the interior node of the current element (left state).
- **$$U_R$$** → Value from the interior node of the neighboring element (right state).

These values form the left and right states needed for a Riemann solver.

**Why This Works:**  
✅ This follows the finite-volume (FV) method's structure, where fluxes are computed at element interfaces.  
✅ Each sub-face naturally forms a FV interface, with two opposing interior states.

---

### **Step 2: Compute the Numerical Flux at the Sub-Face**
Once $$U_L$$ and $$U_R$$ are defined, compute the numerical flux at the sub-face center using any approximate Riemann solver, such as:

$$
F_{\text{sub-face}} = F(U_L, U_R)
$$

**Common Riemann Solvers Used in FV Methods:**
- Roe Flux → Linearized Riemann solver, accurate but requires entropy fix.
- HLL/HLLC Flux → Simple and robust, especially for compressible flows.
- Lax-Friedrichs Flux → Diffusive but stable (can be used in early tests).

🚀 **Result:**  
✅ A first-order accurate FV flux at each sub-face of the global edge.  
✅ Ready for high-order integration using Gauss weights to form a P-th order flux reconstruction at the global edge node.

---

# **2. WENO-Based Flux Reconstruction (Higher-Order)**
For higher-order WENO flux reconstruction, we need to improve the accuracy of $$U_L$$ and $$U_R$$ before computing the numerical flux.

### **Step 1: Form WENO Candidate Reconstructions for $$U_L$$ and $$ U_R $$**
Instead of using just one interior node per element, WENO flux reconstruction uses multiple interior nodes from each element to form a stencil.

For each sub-face of the global edge:
- Define three or more interior nodes from the current element as candidate left states.
- Define three or more interior nodes from the neighboring element as candidate right states.
- Compute smoothness indicators to choose the most accurate state for flux computation.

📌 **How WENO Selects the Best Reconstruction:**
1. Compute multiple candidate reconstructions $$ U_{L,i} $$ and $$ U_{R,i} $$ from different interior nodes.
2. Compute smoothness indicators $$ S_i $$ for each candidate state:

   $$
   S_i = \sum_m \left( \frac{d^m U_i}{dx^m} \right)^2
   $$

3. Compute WENO weights:

   $$
   \alpha_i = \frac{d_i}{(\epsilon + S_i)^p}
   $$

   $$
   \omega_i = \frac{\alpha_i}{\sum_j \alpha_j}
   $$

   where:
   - \( d_i \) are optimal high-order weights.
   - \( \epsilon \) is a small regularization constant.
   - \( p \) controls how aggressively WENO favors smooth stencils.

4. Compute the final reconstructed values:

   $$
   U_L = \sum_i \omega_i U_{L,i}, \quad U_R = \sum_i \omega_i U_{R,i}
   $$

🚀 **Result:**  
✅ Shock-free, high-order accurate flux reconstruction at each sub-face of the global edge.  
✅ Ensures no spurious oscillations near discontinuities.

---

### **Step 2: Compute the Numerical Flux at the Sub-Face**
Now that we have higher-order reconstructed values $$ U_L, U_R $$ from WENO, compute the numerical flux as before:

$$
F_{\text{sub-face}} = F(U_L, U_R)
$$

where \( F(U_L, U_R) \) is evaluated using the same Riemann solver.

🚀 **Result:**  
✅ WENO-reconstructed flux, which is much more accurate than the first-order case.

---

# **3. Final Step: High-Order Edge Flux Reconstruction**
For each **global edge node**, integrate the sub-face fluxes using Gauss quadrature:

$$
F_{\text{edge-node}} = \sum_{i=1}^{P+1} w_i F_{\text{sub-face}, i}
$$

📌 **Why This Works:**  
✅ Recovers full P-th order accuracy at the edge node using high-order integration.  
✅ Maintains robustness near shocks using WENO.  

---

# **Summary of First-Order vs. WENO Reconstruction**
| Approach                               | Accuracy          | Stability           | Computational Cost | Best Used For                   |
| -------------------------------------- | ----------------- | ------------------- | ------------------ | ------------------------------- |
| **First-Order FV Flux Reconstruction** | 1st-order         | Very robust         | Low                | Initial tests, smooth solutions |
| **WENO Flux Reconstruction**           | High-order (P-th) | Handles shocks well | Moderate           | Shock tubes, discontinuities    |

📌 **Your method naturally extends from first-order to high-order WENO, following classical FV techniques.**


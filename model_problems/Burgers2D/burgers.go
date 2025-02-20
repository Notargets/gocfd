package Burgers2D

import (
	"github.com/notargets/gocfd/DG2D"
	"github.com/notargets/gocfd/InputParameters"
	"github.com/notargets/gocfd/utils"
)

/*
The 2D inviscid Burgers' equations can be written in conservative (flux) form as:
∂U∂t+∇⋅F=0
				∂/∂t [ U ] + ∂/∂x [ Fx ] + ∂/∂y [ Fy ] = 0

				U  = [ u ]
                     [ v ]

				Fx = [ ½ u² ] Fy = [  uv  ]
                     [  uv  ]      [ ½ v² ]

Key Properties:
    Hyperbolic system: The fluxes depend nonlinearly on the conserved variables.
    Shock formation: If an initial condition has steep gradients, it will evolve into a shock wave.
    Advection-dominated: It models pure convective transport without viscosity.

At a DG element boundary, we introduce left (L) and right (R) states:
				(u_L, v_L) inside the element
				(u_R, v_R) outside the element

Normal Projection: Reduce to a 1D Riemann Problem

Since flux continuity must be enforced in the normal direction,
we project the velocity onto the normal n=(nx,ny):
				u_n = u * n_x + v * n_y

This defines a 1D projected Burgers' equation along the normal:
				∂u_n/∂t + ∂/∂s (½ u_n²) = 0
where:
				u_n^L = u_L * n_x + v_L * n_y
				u_n^R = u_R * n_x + v_R * n_y

This is now identical to the 1D Burgers’ equation Riemann problem.

2. Solve the 1D Riemann Problem at the Interface

To compute the numerical flux at the interface, we solve the Riemann problem in the normal direction.
Riemann Solution: Shock or Rarefaction

The flux f(u_n) = ½ u_n² satisfies the Rankine-Hugoniot condition with shock
speed s:
				s = (u_n^L + u_n^R) / 2

	Shock Case: (u_n^L > u_n^R), the shock propagates at speed s:
				If s > 0 → use u_n^L
				If s < 0 → use u_n^R
				u_n* =
  					u_n^L,  if s > 0
  					u_n^R,  if s < 0

	Rarefaction Case: (u_n^L < u_n^R), the entropy condition gives:
				u_n* =
  					u_n^L,  if u_n^L > 0
  					0,      if u_n^L < 0 < u_n^R
  					u_n^R,  if u_n^R < 0

3. Compute the Numerical Flux

The numerical flux at the interface is obtained using a Riemann solver, such as Godunov’s exact solver or HLL-type approximate solvers.
Godunov Flux (Exact Riemann Solver)
	Using the exact Riemann solution:
				F*(u_n) = ½ (u_n*)²

	HLL Flux (Approximate Riemann Solver)
				F* = ½ (F_L + F_R) - ½ |s| (u_n^R - u_n^L)
where:
				F_L = ½ (u_n^L)²
				F_R = ½ (u_n^R)²
				s = max(|u_n^L|, |u_n^R|)
    			s is an approximation to the shock speed.

This ensures numerical stability while capturing the correct shock or rarefaction behavior.
*/

type Burgers2D struct {
	U, V utils.Matrix // Dimension Np x K, Np = flux points, inc edges
	// sharing normal flux between elements
	DFR *DG2D.DFR2D
}

func NewBurgers2D(P int) (b2d *Burgers2D) {
	b2d = &Burgers2D{}
	pm := &InputParameters.PlotMeta{}
	b2d.DFR = DG2D.NewDFR2D(P, pm, true)
	return
}

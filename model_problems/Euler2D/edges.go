package Euler2D

import (
	"math"
	"sort"
	"sync"

	"github.com/notargets/gocfd/DG2D"
	"github.com/notargets/gocfd/model_problems/Euler2D/isentropic_vortex"
	"github.com/notargets/gocfd/types"
	"github.com/notargets/gocfd/utils"
)

type EdgeKeySlice []types.EdgeKey

func (p EdgeKeySlice) Len() int           { return len(p) }
func (p EdgeKeySlice) Less(i, j int) bool { return p[i] < p[j] }
func (p EdgeKeySlice) Swap(i, j int)      { p[i], p[j] = p[j], p[i] }

// Sort is a convenience method.
func (p EdgeKeySlice) Sort() { sort.Sort(p) }

func (c *Euler) ParallelSetNormalFluxOnEdges(Time float64) {
	var (
		Ntot = len(c.SortedEdgeKeys)
		wg   = sync.WaitGroup{}
	)
	var ind, end int
	for n := 0; n < c.ParallelDegree; n++ {
		ind, end = c.split1D(Ntot, n)
		wg.Add(1)
		go func(ind, end int) {
			c.SetNormalFluxOnEdges(Time, c.SortedEdgeKeys[ind:end])
			wg.Done()
		}(ind, end)
	}
	wg.Wait()
}

func (c *Euler) SetNormalFluxOnEdges(Time float64, edgeKeys EdgeKeySlice) {
	var (
		dfr                            = c.dfr
		Nint                           = dfr.FluxElement.Nint
		Nedge                          = dfr.FluxElement.Nedge
		Kmax                           = dfr.K
		normalFlux, normalFluxReversed = make([][4]float64, Nedge), make([][4]float64, Nedge)
	)
	for _, en := range edgeKeys {
		e := dfr.Tris.Edges[en]
		//for en, e := range dfr.Tris.Edges {
		switch e.NumConnectedTris {
		case 0:
			panic("unable to handle unconnected edges")
		case 1: // Handle edges with only one triangle - default is edge flux, which will be replaced by a BC flux
			var (
				k           = int(e.ConnectedTris[0])
				edgeNumber  = int(e.ConnectedTriEdgeNumber[0])
				shift       = edgeNumber * Nedge
				riemann     = true
				processFlux bool
				qfD         = Get4DP(c.Q_Face)
			)
			processFlux = true
			normal, _ := c.getEdgeNormal(0, e, en)
			switch e.BCType {
			case types.BC_Far:
				for i := 0; i < Nedge; i++ {
					iL := i + shift
					ind := k + iL*Kmax
					QBC := c.RiemannBC(k, iL, qfD, c.FS.Qinf, normal)
					qfD[0][ind] = QBC[0]
					qfD[1][ind] = QBC[1]
					qfD[2][ind] = QBC[2]
					qfD[3][ind] = QBC[3]
				}
			case types.BC_IVortex:
				// Set the flow variables to the exact solution
				X, Y := c.dfr.FluxX.Data(), c.dfr.FluxY.Data()
				for i := 0; i < Nedge; i++ {
					iL := i + shift
					indFlux := k + (2*Nint+iL)*Kmax
					x, y := X[indFlux], Y[indFlux]
					iv := c.AnalyticSolution.(*isentropic_vortex.IVortex)
					rho, rhoU, rhoV, E := iv.GetStateC(Time, x, y)
					ind := k + iL*Kmax
					var QBC [4]float64
					if riemann {
						Qinf := [4]float64{rho, rhoU, rhoV, E}
						QBC = c.RiemannBC(k, iL, qfD, Qinf, normal)
					} else {
						QBC = [4]float64{rho, rhoU, rhoV, E}
					}
					qfD[0][ind] = QBC[0]
					qfD[1][ind] = QBC[1]
					qfD[2][ind] = QBC[2]
					qfD[3][ind] = QBC[3]
				}
			case types.BC_Wall, types.BC_Cyl:
				processFlux = false
				for i := 0; i < Nedge; i++ {
					ie := i + shift
					p := c.FS.GetFlowFunctionAtIndex(k+ie*Kmax, qfD, StaticPressure)
					for n := 0; n < 4; n++ {
						normalFlux[i][0] = 0
						normalFlux[i][3] = 0
						normalFlux[i][1] = normal[0] * p
						normalFlux[i][2] = normal[1] * p
					}
				}
				c.SetNormalFluxOnRTEdge(k, edgeNumber, normalFlux, e.IInII[0])
			case types.BC_PeriodicReversed, types.BC_Periodic:
				processFlux = false
			}
			if processFlux {
				var Fx, Fy [4]float64
				for i := 0; i < Nedge; i++ {
					ie := i + shift
					Fx, Fy = c.CalculateFlux(k, ie, qfD)
					for n := 0; n < 4; n++ {
						normalFlux[i][n] = normal[0]*Fx[n] + normal[1]*Fy[n]
					}
				}
				c.SetNormalFluxOnRTEdge(k, edgeNumber, normalFlux, e.IInII[0])
			}
		case 2: // Handle edges with two connected tris - shared faces
			var (
				kL, kR                   = int(e.ConnectedTris[0]), int(e.ConnectedTris[1])
				edgeNumberL, edgeNumberR = int(e.ConnectedTriEdgeNumber[0]), int(e.ConnectedTriEdgeNumber[1])
				shiftL, shiftR           = edgeNumberL * Nedge, edgeNumberR * Nedge
				fluxLeft, fluxRight      [2][4]float64
				qfD                      = Get4DP(c.Q_Face)
			)
			switch c.FluxCalcAlgo {
			case FLUX_Average:
				averageFluxN := func(f1, f2 [2][4]float64, normal [2]float64) (fnorm [4]float64, fnormR [4]float64) {
					var (
						fave [2][4]float64
					)
					for n := 0; n < 4; n++ {
						for ii := 0; ii < 2; ii++ {
							fave[ii][n] = 0.5 * (f1[ii][n] + f2[ii][n])
						}
						fnorm[n] = normal[0]*fave[0][n] + normal[1]*fave[1][n]
						fnormR[n] = -fnorm[n]
					}
					return
				}
				normal, _ := c.getEdgeNormal(0, e, en)
				for i := 0; i < Nedge; i++ {
					iL := i + shiftL
					iR := Nedge - 1 - i + shiftR // Shared edges run in reverse order relative to each other
					fluxLeft[0], fluxLeft[1] = c.CalculateFlux(kL, iL, qfD)
					fluxRight[0], fluxRight[1] = c.CalculateFlux(kR, iR, qfD) // Reverse the right edge to match
					normalFlux[i], normalFluxReversed[Nedge-1-i] = averageFluxN(fluxLeft, fluxRight, normal)
				}
			case FLUX_LaxFriedrichs:
				var (
					rhoL, uL, vL, pL, CL float64
					rhoR, uR, vR, pR, CR float64
					maxV                 float64
				)
				maxVF := func(u, v, p, rho, C float64) (vmax float64) {
					vmax = math.Sqrt(u*u+v*v) + C
					return
				}
				normal, _ := c.getEdgeNormal(0, e, en)
				for i := 0; i < Nedge; i++ {
					iL := i + shiftL
					iR := Nedge - 1 - i + shiftR // Shared edges run in reverse order relative to each other
					qL := c.GetQQ(kL, iL, qfD)
					rhoL, uL, vL = qL[0], qL[1]/qL[0], qL[2]/qL[0]
					pL, CL = c.FS.GetFlowFunction(qL, StaticPressure), c.FS.GetFlowFunction(qL, SoundSpeed)
					qR := c.GetQQ(kR, iR, qfD)
					rhoR, uR, vR = qR[0], qR[1]/qR[0], qR[2]/qR[0]
					pR, CR = c.FS.GetFlowFunction(qR, StaticPressure), c.FS.GetFlowFunction(qR, SoundSpeed)
					fluxLeft[0], fluxLeft[1] = c.CalculateFlux(kL, iL, qfD)
					fluxRight[0], fluxRight[1] = c.CalculateFlux(kR, iR, qfD) // Reverse the right edge to match
					maxV = math.Max(maxVF(uL, vL, pL, rhoL, CL), maxVF(uR, vR, pR, rhoR, CR))
					indL, indR := kL+iL*Kmax, kR+iR*Kmax
					for n := 0; n < 4; n++ {
						normalFlux[i][n] = 0.5 * (normal[0]*(fluxLeft[0][n]+fluxRight[0][n]) + normal[1]*(fluxLeft[1][n]+fluxRight[1][n]))
						normalFlux[i][n] += 0.5 * maxV * (qfD[n][indL] - qfD[n][indR])
						normalFluxReversed[Nedge-1-i][n] = -normalFlux[i][n]
					}
				}
			case FLUX_Roe:
				var (
					rhoL, uL, vL, pL float64
					rhoR, uR, vR, pR float64
					hL, hR           float64
					Gamma            = c.FS.Gamma
					GM1              = Gamma - 1
					qfD              = Get4DP(c.Q_Face)
				)
				normal, _ := c.getEdgeNormal(0, e, en)
				rotateMomentum := func(k, i int) {
					ind := k + i*Kmax
					um, vm := qfD[1][ind], qfD[2][ind]
					qfD[1][ind] = um*normal[0] + vm*normal[1]
					qfD[2][ind] = -um*normal[1] + vm*normal[0]
				}
				for i := 0; i < Nedge; i++ {
					iL := i + shiftL
					iR := Nedge - 1 - i + shiftR // Shared edges run in reverse order relative to each other
					// Rotate the momentum into face normal coordinates before calculating fluxes
					rotateMomentum(kL, iL)
					rotateMomentum(kR, iR)
					qL := c.GetQQ(kL, iL, qfD)
					rhoL, uL, vL = qL[0], qL[1]/qL[0], qL[2]/qL[0]
					pL = c.FS.GetFlowFunction(qL, StaticPressure)
					qR := c.GetQQ(kR, iR, qfD)
					rhoR, uR, vR = qR[0], qR[1]/qR[0], qR[2]/qR[0]
					pR = c.FS.GetFlowFunction(qR, StaticPressure)
					fluxLeft[0], _ = c.CalculateFlux(kL, iL, qfD)
					fluxRight[0], _ = c.CalculateFlux(kR, iR, qfD) // Reverse the right edge to match
					/*
					   HM = (EnerM+pM).dd(rhoM);  HP = (EnerP+pP).dd(rhoP);
					*/
					// Enthalpy
					hL, hR = c.FS.GetFlowFunction(qL, Enthalpy), c.FS.GetFlowFunction(qR, Enthalpy)
					/*
						rhoMs = sqrt(rhoM); rhoPs = sqrt(rhoP);
						rhoMsPs = rhoMs + rhoPs;

						rho = rhoMs.dm(rhoPs);
						u   = (rhoMs.dm(uM) + rhoPs.dm(uP)).dd(rhoMsPs);
						v   = (rhoMs.dm(vM) + rhoPs.dm(vP)).dd(rhoMsPs);
						H   = (rhoMs.dm(HM) + rhoPs.dm(HP)).dd(rhoMsPs);
						c2 = gm1 * (H - 0.5*(sqr(u)+sqr(v)));
						c = sqrt(c2);
					*/
					// Compute Roe average variables
					rhoLs, rhoRs := math.Sqrt(rhoL), math.Sqrt(rhoR)
					rhoLsRs := rhoLs + rhoRs

					rho := rhoLs * rhoRs
					u := (rhoLs*uL + rhoRs*uR) / rhoLsRs
					v := (rhoLs*vL + rhoRs*vR) / rhoLsRs
					h := (rhoLs*hL + rhoRs*hR) / rhoLsRs
					c2 := GM1 * (h - 0.5*(u*u+v*v))
					c := math.Sqrt(c2)
					/*
					   dW1 = -0.5*(rho.dm(uP-uM)).dd(c) + 0.5*(pP-pM).dd(c2);
					   dW2 = (rhoP-rhoM) - (pP-pM).dd(c2);
					   dW3 = rho.dm(vP-vM);
					   dW4 = 0.5*(rho.dm(uP-uM)).dd(c) + 0.5*(pP-pM).dd(c2);

					   dW1 = abs(u-c).dm(dW1);
					   dW2 = abs(u  ).dm(dW2);
					   dW3 = abs(u  ).dm(dW3);
					   dW4 = abs(u+c).dm(dW4);
					*/
					// Riemann fluxes
					dW1 := -0.5*(rho*(uR-uL))/c + 0.5*(pR-pL)/c2
					dW2 := (rhoR - rhoL) - (pR-pL)/c2
					dW3 := rho * (vR - vL)
					dW4 := 0.5*(rho*(uR-uL))/c + 0.5*(pR-pL)/c2
					dW1 = math.Abs(u-c) * dW1
					dW2 = math.Abs(u) * dW2
					dW3 = math.Abs(u) * dW3
					dW4 = math.Abs(u+c) * dW4
					/*
					   DMat fx = (fxQP+fxQM)/2.0;
					   fx(All,1) -= (dW1               + dW2                                   + dW4              )/2.0;
					   fx(All,2) -= (dW1.dm(u-c)       + dW2.dm(u)                             + dW4.dm(u+c)      )/2.0;
					   fx(All,3) -= (dW1.dm(v)         + dW2.dm(v)                 + dW3       + dW4.dm(v)        )/2.0;
					   fx(All,4) -= (dW1.dm(H-u.dm(c)) + dW2.dm(sqr(u)+sqr(v))/2.0 + dW3.dm(v) + dW4.dm(H+u.dm(c)))/2.0;
					*/
					// Form Roe Fluxes
					for n := 0; n < 4; n++ {
						normalFlux[i][n] = 0.5 * (fluxLeft[0][n] + fluxRight[0][n]) // Ave of normal component of flux
					}
					normalFlux[i][0] -= 0.5 * (dW1 + dW2 + dW4)
					normalFlux[i][1] -= 0.5 * (dW1*(u-c) + dW2*u + dW4*(u+c))
					normalFlux[i][2] -= 0.5 * (dW1*v + dW2*v + dW3 + dW4*v)
					normalFlux[i][3] -= 0.5 * (dW1*(h-u*c) + 0.5*dW2*(u*u+v*v) + dW3*v + dW4*(h+u*c))
					/*
					   flux = fx;    fx2.borrow(Ngf, fx.pCol(2)); fx3.borrow(Ngf, fx.pCol(3));
					   flux(All,2) = lnx.dm(fx2) - lny.dm(fx3);
					   flux(All,3) = lny.dm(fx2) + lnx.dm(fx3);
					*/
					// rotate back to Cartesian
					normalFlux[i][1], normalFlux[i][2] = normal[0]*normalFlux[i][1]-normal[1]*normalFlux[i][2], normal[1]*normalFlux[i][1]+normal[0]*normalFlux[i][2]
					for n := 0; n < 4; n++ {
						normalFluxReversed[Nedge-1-i][n] = -normalFlux[i][n]
					}
				}
			}
			c.SetNormalFluxOnRTEdge(kL, edgeNumberL, normalFlux, e.IInII[0])
			c.SetNormalFluxOnRTEdge(kR, edgeNumberR, normalFluxReversed, e.IInII[1])
		}
	}
	return
}

func (c *Euler) getEdgeNormal(conn int, e *DG2D.Edge, en types.EdgeKey) (normal, scaledNormal [2]float64) {
	var (
		dfr = c.dfr
	)
	norm := func(vec [2]float64) (n float64) {
		n = math.Sqrt(vec[0]*vec[0] + vec[1]*vec[1])
		return
	}
	normalize := func(vec [2]float64) (normed [2]float64) {
		n := norm(vec)
		for i := 0; i < 2; i++ {
			normed[i] = vec[i] / n
		}
		return
	}
	revDir := bool(e.ConnectedTriDirection[conn])
	x1, x2 := DG2D.GetEdgeCoordinates(en, revDir, dfr.VX, dfr.VY)
	dx := [2]float64{x2[0] - x1[0], x2[1] - x1[1]}
	normal = normalize([2]float64{-dx[1], dx[0]})
	scaledNormal[0] = normal[0] * e.IInII[conn]
	scaledNormal[1] = normal[1] * e.IInII[conn]
	return
}

func (c *Euler) ProjectFluxToEdge(edgeFlux [][2][4]float64, e *DG2D.Edge, en types.EdgeKey) {
	/*
				Projects a 2D flux: [F, G] onto the face normal, then multiplies by the element/edge rqtio of normals, ||n||
		 		And places the scaled and projected normal flux into the degree of freedom F_RT_DOF
	*/
	var (
		dfr        = c.dfr
		Nedge      = dfr.FluxElement.Nedge
		Nint       = dfr.FluxElement.Nint
		Kmax       = dfr.K
		conn       = 0
		k          = int(e.ConnectedTris[conn])
		edgeNumber = int(e.ConnectedTriEdgeNumber[conn])
		shift      = edgeNumber * Nedge
	)
	// Get scaling factor ||n|| for each edge, multiplied by untransformed normals
	_, scaledNormal := c.getEdgeNormal(conn, e, en)
	for i := 0; i < Nedge; i++ {
		iL := i + shift
		// Project the flux onto the scaled scaledNormal
		for n := 0; n < 4; n++ {
			// Place normed/scaled flux into the RT element space
			rtD := c.F_RT_DOF[n].Data()
			ind := k + (2*Nint+iL)*Kmax
			rtD[ind] = scaledNormal[0]*edgeFlux[i][0][n] + scaledNormal[1]*edgeFlux[i][1][n]
		}
	}
}

func (c *Euler) SetNormalFluxOnRTEdge(k, edgeNumber int, edgeNormalFlux [][4]float64, IInII float64) {
	/*
		Takes the normal flux (aka "projected flux") multiplies by the ||n|| ratio of edge normals and sets that value for
		the F_RT_DOF degree of freedom locations for this [k, edgenumber] group
	*/
	var (
		dfr   = c.dfr
		Nedge = dfr.FluxElement.Nedge
		Nint  = dfr.FluxElement.Nint
		Kmax  = dfr.K
		shift = edgeNumber * Nedge
	)
	// Get scaling factor ||n|| for each edge, multiplied by untransformed normals
	for n := 0; n < 4; n++ {
		rtD := c.F_RT_DOF[n].Data()
		for i := 0; i < Nedge; i++ {
			// Place normed/scaled flux into the RT element space
			ind := k + (2*Nint+i+shift)*Kmax
			rtD[ind] = edgeNormalFlux[i][n] * IInII
		}
	}
}

func (c *Euler) InterpolateSolutionToEdges(Q [4]utils.Matrix) {
	// Interpolate from solution points to edges using precomputed interpolation matrix
	var wg = sync.WaitGroup{}
	for n := 0; n < 4; n++ {
		wg.Add(1)
		go func(n int) {
			c.Q_Face[n] = c.dfr.FluxEdgeInterpMatrix.Mul(Q[n])
			wg.Done()
		}(n)
	}
	wg.Wait()
}

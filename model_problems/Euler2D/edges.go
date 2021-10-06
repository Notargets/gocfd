package Euler2D

import (
	"math"
	"sort"

	"github.com/notargets/gocfd/types"
	"github.com/notargets/gocfd/utils"
)

type EdgeKeySlice []types.EdgeKey

func (p EdgeKeySlice) Len() int { return len(p) }
func (p EdgeKeySlice) Less(i, j int) bool {
	// Sorted to produce groups of vertex centered edges on the right edge vertex
	vnode1, vnode2 := int(p[i]>>32), int(p[j]>>32)
	return vnode1 < vnode2
}
func (p EdgeKeySlice) Swap(i, j int) { p[i], p[j] = p[j], p[i] }

// Sort is a convenience method.
func (p EdgeKeySlice) Sort() { sort.Sort(p) }

type EdgeKeySliceSortLeft EdgeKeySlice

func (p EdgeKeySliceSortLeft) Sort() { sort.Sort(p) }

func (p EdgeKeySliceSortLeft) Len() int      { return len(p) }
func (p EdgeKeySliceSortLeft) Swap(i, j int) { p[i], p[j] = p[j], p[i] }
func (p EdgeKeySliceSortLeft) Less(i, j int) bool {
	getLeftVert := func(ek types.EdgeKey) (left int) {
		enTmp := ek >> 32
		left = int(ek - enTmp*(1<<32))
		return
	}
	// Sorted to produce groups of vertex centered edges on the left edge vertex
	vnode1, vnode2 := getLeftVert(p[i]), getLeftVert(p[j])
	return vnode1 < vnode2
}

func (c *Euler) SetNormalFluxOnEdges(Time float64, CalculateDT bool, Jdet, DT []utils.Matrix, F_RT_DOF, Q_Face [][4]utils.Matrix, edgeKeys EdgeKeySlice, EdgeQ1, EdgeQ2 [][4]float64) (waveSpeedMax float64) {
	var (
		Nedge                          = c.dfr.FluxElement.Nedge
		normalFlux, normalFluxReversed = EdgeQ1, EdgeQ2
		pm                             = c.Partitions
	)
	for _, en := range edgeKeys {
		e := c.dfr.Tris.Edges[en]
		switch e.NumConnectedTris {
		case 0:
			panic("unable to handle unconnected edges")
		case 1: // Handle edges with only one triangle - default is edge flux, which will be replaced by a BC flux
			var (
				k, Kmax, bn         = pm.GetLocalK(int(e.ConnectedTris[0]))
				edgeNumber          = int(e.ConnectedTriEdgeNumber[0])
				shift               = edgeNumber * Nedge
				calculateNormalFlux bool
			)
			calculateNormalFlux = true
			//normal, _ := c.getEdgeNormal(0, e, en)
			normal := e.GetEdgeNormal(0, en, c.dfr)
			switch e.BCType {
			case types.BC_Far:
				c.FarBC(k, Kmax, shift, Q_Face[bn], normal)
			case types.BC_IVortex:
				c.IVortexBC(Time, k, Kmax, shift, Q_Face[bn], normal)
			case types.BC_Wall, types.BC_Cyl:
				calculateNormalFlux = false
				c.WallBC(k, Kmax, Q_Face[bn], shift, normal, normalFlux) // Calculates normal flux directly
			case types.BC_PeriodicReversed, types.BC_Periodic:
				// One edge of the Periodic BC leads to calculation of both sides within the connected tris section, so noop here
				continue
			}
			if calculateNormalFlux {
				var Fx, Fy [4]float64
				for i := 0; i < Nedge; i++ {
					ie := i + shift
					ind := k + ie*Kmax
					Fx, Fy = c.CalculateFlux(Q_Face[bn], ind)
					for n := 0; n < 4; n++ {
						normalFlux[i][n] = normal[0]*Fx[n] + normal[1]*Fy[n]
					}
				}
			}
			c.SetNormalFluxOnRTEdge(k, Kmax, F_RT_DOF[bn], edgeNumber, normalFlux, e.IInII[0])
		case 2: // Handle edges with two connected tris - shared faces
			var (
				//				kL, kR                   = int(e.ConnectedTris[0]), int(e.ConnectedTris[1])
				kL, KmaxL, bnL           = pm.GetLocalK(int(e.ConnectedTris[0]))
				kR, KmaxR, bnR           = pm.GetLocalK(int(e.ConnectedTris[1]))
				edgeNumberL, edgeNumberR = int(e.ConnectedTriEdgeNumber[0]), int(e.ConnectedTriEdgeNumber[1])
				shiftL, shiftR           = edgeNumberL * Nedge, edgeNumberR * Nedge
			)
			//normal, _ := c.getEdgeNormal(0, e, en)
			normal := e.GetEdgeNormal(0, en, c.dfr)
			switch c.FluxCalcAlgo {
			case FLUX_Average:
				c.AvgFlux(kL, kR, KmaxL, KmaxR, shiftL, shiftR,
					Q_Face[bnL], Q_Face[bnR], normal, normalFlux, normalFluxReversed)
			case FLUX_LaxFriedrichs:
				c.LaxFlux(kL, kR, KmaxL, KmaxR, shiftL, shiftR,
					Q_Face[bnL], Q_Face[bnR], normal, normalFlux, normalFluxReversed)
			case FLUX_Roe:
				c.RoeFlux(kL, kR, KmaxL, KmaxR, shiftL, shiftR,
					Q_Face[bnL], Q_Face[bnR], normal, normalFlux, normalFluxReversed)
			case FLUX_RoeER:
				c.RoeERFlux(kL, kR, KmaxL, KmaxR, shiftL, shiftR,
					Q_Face[bnL], Q_Face[bnR], normal, normalFlux, normalFluxReversed)
			}
			c.SetNormalFluxOnRTEdge(kL, KmaxL, F_RT_DOF[bnL], edgeNumberL, normalFlux, e.IInII[0])
			c.SetNormalFluxOnRTEdge(kR, KmaxR, F_RT_DOF[bnR], edgeNumberR, normalFluxReversed, e.IInII[1])
		}
		if CalculateDT {
			var (
				Np1  = c.dfr.N + 1
				Np12 = float64(Np1 * Np1)
			)
			var (
				conn        = 0
				edgeNum     = int(e.ConnectedTriEdgeNumber[conn])
				shift       = edgeNum * Nedge
				edgeLen     = e.GetEdgeLength()
				k, Kmax, bn = pm.GetLocalK(int(e.ConnectedTris[conn]))
			)
			// fmt.Printf("N, Np12, edgelen, Jdet = %d,%8.5f,%8.5f,%8.5f\n", c.dfr.N, Np12, edgeLen, Jdet)
			fs := 0.5 * Np12 * edgeLen / Jdet[bn].DataP[k]
			edgeMax := -100.
			for i := shift; i < shift+Nedge; i++ {
				ind := k + i*Kmax
				C := c.FS.GetFlowFunction(Q_Face[bn], ind, SoundSpeed)
				U := c.FS.GetFlowFunction(Q_Face[bn], ind, Velocity)
				waveSpeed := fs * (U + C)
				waveSpeedMax = math.Max(waveSpeed, waveSpeedMax)
				if waveSpeed > edgeMax {
					edgeMax = waveSpeed
				}
			}
			if edgeMax > DT[bn].DataP[k] {
				DT[bn].DataP[k] = edgeMax
			}
			if e.NumConnectedTris == 2 { // Add the wavespeed to the other tri connected to this edge if needed
				k, Kmax, bn = pm.GetLocalK(int(e.ConnectedTris[1]))
				if edgeMax > DT[bn].DataP[k] {
					DT[bn].DataP[k] = edgeMax
				}
			}
		}
	}
	return
}

func (c *Euler) SetNormalFluxOnRTEdge(k, Kmax int, F_RT_DOF [4]utils.Matrix, edgeNumber int, edgeNormalFlux [][4]float64, IInII float64) {
	/*
		Takes the normal flux (aka "projected flux") multiplies by the ||n|| ratio of edge normals and sets that value for
		the F_RT_DOF degree of freedom locations for this [k, edgenumber] group
	*/
	var (
		dfr   = c.dfr
		Nedge = dfr.FluxElement.Nedge
		Nint  = dfr.FluxElement.Nint
		shift = edgeNumber * Nedge
	)
	// Get scaling factor ||n|| for each edge, multiplied by untransformed normals
	for n := 0; n < 4; n++ {
		rtD := F_RT_DOF[n].DataP
		for i := 0; i < Nedge; i++ {
			// Place normed/scaled flux into the RT element space
			ind := k + (2*Nint+i+shift)*Kmax
			rtD[ind] = edgeNormalFlux[i][n] * IInII
		}
	}
}

func (c *Euler) InterpolateSolutionToEdges(Q, Q_Face [4]utils.Matrix) {
	// Interpolate from solution points to edges using precomputed interpolation matrix
	for n := 0; n < 4; n++ {
		c.dfr.FluxEdgeInterp.Mul(Q[n], Q_Face[n])
	}
	return
}

package Euler2D

import (
	"math"
	"sort"

	"github.com/notargets/gocfd/DG2D"

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

type EdgeValueStorage struct {
	EdgeFluxStorage     [4]utils.Matrix       // Normal flux storage, dimension is NpEdge x Nedges
	EdgeSolutionStorage [4]utils.Matrix       // Edge Solution storage, dimension is NpEdge x Nedges
	StorageIndex        map[types.EdgeKey]int // Index into normal flux storage using edge key
}

type ValueType uint8

const (
	EdgeFlux ValueType = iota
	SolutionValues
)

func (c *Euler) NewEdgeStorage() (nf *EdgeValueStorage) {
	nf = &EdgeValueStorage{
		StorageIndex: make(map[types.EdgeKey]int),
	}
	// Allocate memory for Normal flux
	NumEdges := len(c.dfr.Tris.Edges)
	for n := 0; n < 4; n++ {
		nf.EdgeFluxStorage[n] = utils.NewMatrix(NumEdges, c.dfr.FluxElement.Nedge)
		nf.EdgeSolutionStorage[n] = utils.NewMatrix(NumEdges, c.dfr.FluxElement.Nedge)
	}
	var index int
	for en, _ := range c.dfr.Tris.Edges {
		nf.StorageIndex[en] = index
		index++
	}
	return
}

func (nf *EdgeValueStorage) GetEdgeValues(valType ValueType, kGlobal, localEdgeNumber int, dfr *DG2D.DFR2D) (EdgeValues [4][]float64, sign int) {
	var (
		k       = kGlobal
		Kmax    = dfr.K
		edgeNum = localEdgeNumber
		Nedge   = dfr.FluxElement.Nedge
		target  [4]utils.Matrix
	)
	switch valType {
	case EdgeFlux:
		target = nf.EdgeFluxStorage
	case SolutionValues:
		target = nf.EdgeSolutionStorage
	}
	ind := k + Kmax*edgeNum
	en := dfr.EdgeNumber[ind]
	edgeIndex := nf.StorageIndex[en]
	ind = edgeIndex * Nedge
	for n := 0; n < 4; n++ {
		EdgeValues[n] = target[n].DataP[ind : ind+Nedge]
	}
	if int(dfr.Tris.Edges[en].ConnectedTris[0]) == kGlobal {
		// These values were stored in this element's order
		sign = 1
	} else {
		// These values should be reversed in direction (and if normals, sign as well)
		sign = -1
	}
	return
}

func (c *Euler) CalculateNormalFlux(Time float64, CalculateDT bool, Jdet, DT []utils.Matrix, Q_Face [][4]utils.Matrix,
	edgeKeys EdgeKeySlice, EdgeQ1, EdgeQ2 [][4]float64) (waveSpeedMax float64) {
	var (
		Nedge       = c.dfr.FluxElement.Nedge
		normalFlux  = EdgeQ1
		avgSolution = EdgeQ2
		pm          = c.Partitions
		KmaxGlobal  = c.dfr.K
	)
	for _, en := range edgeKeys {
		e := c.dfr.Tris.Edges[en]
		var (
			k0Global    = int(e.ConnectedTris[0])
			edgeNumber0 = int(e.ConnectedTriEdgeNumber[0])
			faceInd     = k0Global + KmaxGlobal*edgeNumber0
			normal0     = [2]float64{c.dfr.FaceNorm[0].DataP[faceInd], c.dfr.FaceNorm[1].DataP[faceInd]}
		)
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
			switch e.BCType {
			case types.BC_Far:
				c.FarBC(k, Kmax, shift, Q_Face[bn], normal0)
			case types.BC_IVortex:
				c.IVortexBC(Time, k, Kmax, shift, Q_Face[bn], normal0)
			case types.BC_Wall, types.BC_Cyl:
				calculateNormalFlux = false
				c.WallBC(k, Kmax, Q_Face[bn], shift, normal0, normalFlux) // Calculates normal flux directly
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
						normalFlux[i][n] = normal0[0]*Fx[n] + normal0[1]*Fy[n]
					}
				}
			}
			// Store the average solution for this edge - there's only one tri, so it's just a copy of this edge
			for i := 0; i < Nedge; i++ {
				ie := i + shift
				ind := k + ie*Kmax
				for n := 0; n < 4; n++ {
					avgSolution[i][n] = Q_Face[bn][n].DataP[ind]
				}
			}
		case 2: // Handle edges with two connected tris - shared faces
			var (
				//				kL, kR                   = int(e.ConnectedTris[0]), int(e.ConnectedTris[1])
				kL, KmaxL, bnL           = pm.GetLocalK(int(e.ConnectedTris[0]))
				kR, KmaxR, bnR           = pm.GetLocalK(int(e.ConnectedTris[1]))
				edgeNumberL, edgeNumberR = int(e.ConnectedTriEdgeNumber[0]), int(e.ConnectedTriEdgeNumber[1])
				shiftL, shiftR           = edgeNumberL * Nedge, edgeNumberR * Nedge
			)
			switch c.FluxCalcAlgo {
			case FLUX_Average:
				c.AvgFlux(kL, kR, KmaxL, KmaxR, shiftL, shiftR, Q_Face[bnL], Q_Face[bnR], normal0, normalFlux)
			case FLUX_LaxFriedrichs:
				c.LaxFlux(kL, kR, KmaxL, KmaxR, shiftL, shiftR, Q_Face[bnL], Q_Face[bnR], normal0, normalFlux)
			case FLUX_Roe:
				c.RoeFlux(kL, kR, KmaxL, KmaxR, shiftL, shiftR, Q_Face[bnL], Q_Face[bnR], normal0, normalFlux)
			case FLUX_RoeER:
				c.RoeERFlux(kL, kR, KmaxL, KmaxR, shiftL, shiftR, Q_Face[bnL], Q_Face[bnR], normal0, normalFlux)
			}
			// Store the average solution for this edge - average of the two connected tri edges
			for i := 0; i < Nedge; i++ {
				indL := kL + (i+shiftL)*KmaxL
				indR := kR + (Nedge-1-i+shiftR)*KmaxR
				for n := 0; n < 4; n++ {
					avgSolution[i][n] = 0.5 * (Q_Face[bnL][n].DataP[indL] + Q_Face[bnR][n].DataP[indR])
				}
			}
		}
		// Load the normal flux into the global normal flux storage
		edgeIndex := c.EdgeStore.StorageIndex[en]
		for i := 0; i < Nedge; i++ {
			ind := i + edgeIndex*Nedge
			for n := 0; n < 4; n++ {
				c.EdgeStore.EdgeFluxStorage[n].DataP[ind] = normalFlux[i][n]
				c.EdgeStore.EdgeSolutionStorage[n].DataP[ind] = avgSolution[i][n]
			}
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

func (c *Euler) SetRTFluxOnEdges(myThread, Kmax int, F_RT_DOF [4]utils.Matrix) {
	var (
		dfr        = c.dfr
		Nedge      = dfr.FluxElement.Nedge
		Nint       = dfr.FluxElement.Nint
		KmaxGlobal = c.dfr.K
	)
	for k := 0; k < Kmax; k++ {
		kGlobal := c.Partitions.GetGlobalK(k, myThread)
		for edgeNum := 0; edgeNum < 3; edgeNum++ {
			shift := edgeNum * Nedge
			//nFlux, sign := c.EdgeStore.GetEdgeNormalFlux(kGlobal, edgeNum, dfr)
			nFlux, sign := c.EdgeStore.GetEdgeValues(EdgeFlux, kGlobal, edgeNum, dfr)
			ind2 := kGlobal + KmaxGlobal*edgeNum
			IInII := dfr.IInII.DataP[ind2]
			for n := 0; n < 4; n++ {
				rtD := F_RT_DOF[n].DataP
				for i := 0; i < Nedge; i++ {
					// Place normed/scaled flux into the RT element space
					ind := k + (2*Nint+i+shift)*Kmax
					if sign > 0 {
						rtD[ind] = nFlux[n][i] * IInII
					} else {
						rtD[ind] = -nFlux[n][Nedge-i-1] * IInII
					}
				}
			}
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

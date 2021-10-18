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
	EdgeFluxStorage                      [4]utils.Matrix       // Normal flux storage, dimension is NpEdge x Nedges
	EdgeSolutionStorage                  [4]utils.Matrix       // Edge Solution storage, dimension is NpEdge x Nedges
	EdgeArtificialDissipationFluxStorage [4]utils.Matrix       // Stores the numerical flux of gradient of U(X,Y)
	StorageIndex                         map[types.EdgeKey]int // Index into normal flux storage using edge key
	PMap                                 *PartitionMap
}

type ValueType uint8

const (
	EdgeFlux ValueType = iota
	SolutionValues
	ArtificialDissipationFlux
)

func (c *Euler) NewEdgeStorage() (nf *EdgeValueStorage) {
	nf = &EdgeValueStorage{
		StorageIndex: make(map[types.EdgeKey]int),
		PMap:         c.Partitions,
	}
	// Allocate memory for Normal flux
	NumEdges := len(c.dfr.Tris.Edges)
	for n := 0; n < 4; n++ {
		nf.EdgeFluxStorage[n] = utils.NewMatrix(NumEdges, c.dfr.FluxElement.Nedge)
		nf.EdgeSolutionStorage[n] = utils.NewMatrix(NumEdges, c.dfr.FluxElement.Nedge)
		nf.EdgeArtificialDissipationFluxStorage[n] = utils.NewMatrix(NumEdges, c.dfr.FluxElement.Nedge)
	}
	var index int
	for np := 0; np < c.Partitions.ParallelDegree; np++ {
		for _, en := range c.SortedEdgeKeys[np] {
			nf.StorageIndex[en] = index
			index++
		}
	}
	return
}

func (nf *EdgeValueStorage) GetEdgeValues(valType ValueType, myThread, kLocal, varNum, localEdgeNumber int, dfr *DG2D.DFR2D) (EdgeValues []float64, sign int) {
	var (
		kGlobal = nf.PMap.GetGlobalK(kLocal, myThread)
		Kmax    = dfr.K
		edgeNum = localEdgeNumber
		Nedge   = dfr.FluxElement.Nedge
		target  utils.Matrix
	)
	switch valType {
	case EdgeFlux:
		target = nf.EdgeFluxStorage[varNum]
	case SolutionValues:
		target = nf.EdgeSolutionStorage[varNum]
	case ArtificialDissipationFlux:
		target = nf.EdgeArtificialDissipationFluxStorage[varNum]
	}
	ind := kGlobal + Kmax*edgeNum
	en := dfr.EdgeNumber[ind]
	edgeIndex := nf.StorageIndex[en]
	ind = edgeIndex * Nedge
	EdgeValues = target.DataP[ind : ind+Nedge]
	if int(dfr.Tris.Edges[en].ConnectedTris[0]) == kGlobal {
		// These values were stored in this element's order
		sign = 1
	} else {
		// These values should be reversed in direction (and if normals, sign as well)
		sign = -1
	}
	return
}

func (c *Euler) GetFaceNormal(kGlobal, edgeNumber int) (normal [2]float64) {
	var (
		KmaxGlobal = c.dfr.K
		faceInd    = kGlobal + KmaxGlobal*edgeNumber
	)
	normal = [2]float64{c.dfr.FaceNorm[0].DataP[faceInd], c.dfr.FaceNorm[1].DataP[faceInd]}
	return
}

func (c *Euler) CalculateNormalFlux(Time float64, CalculateDT bool, Jdet, DT []utils.Matrix, Q_Face [][4]utils.Matrix,
	edgeKeys EdgeKeySlice, EdgeQ1, EdgeQ2 [][4]float64) (waveSpeedMax float64) {
	var (
		Nedge       = c.dfr.FluxElement.Nedge
		normalFlux  = EdgeQ1
		avgSolution = EdgeQ2
		pm          = c.Partitions
	)
	for _, en := range edgeKeys {
		e := c.dfr.Tris.Edges[en]
		var (
			k0Global             = int(e.ConnectedTris[0])
			k0, Kmax0, myThread0 = pm.GetLocalK(int(e.ConnectedTris[0]))
			edgeNumber0          = int(e.ConnectedTriEdgeNumber[0])
			normal0              = c.GetFaceNormal(k0Global, edgeNumber0)
		)
		switch e.NumConnectedTris {
		case 0:
			panic("unable to handle unconnected edges")
		case 1: // Handle edges with only one triangle - default is edge flux, which will be replaced by a BC flux
			c.calculateBCs(e, Nedge, Time,
				k0, Kmax0, edgeNumber0, myThread0,
				normal0, normalFlux, avgSolution, Q_Face)
		case 2: // Handle edges with two connected tris - shared faces
			var (
				kR, KmaxR, myThreadR = pm.GetLocalK(int(e.ConnectedTris[1]))
				edgeNumberR          = int(e.ConnectedTriEdgeNumber[1])
			)
			c.calculateNumericalFlux(Nedge,
				k0, Kmax0, edgeNumber0, myThread0,
				kR, KmaxR, edgeNumberR, myThreadR,
				normal0, normalFlux, avgSolution, Q_Face)
		}
		if CalculateDT {
			waveSpeedMax = c.calculateLocalDT(e, Nedge, Q_Face, Jdet, DT)
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
	}
	return
}

func (c *Euler) calculateLocalDT(e *DG2D.Edge, Nedge int,
	Q_Face [][4]utils.Matrix, Jdet, DT []utils.Matrix) (waveSpeedMax float64) {
	var (
		pm                = c.Partitions
		edgeNum           = int(e.ConnectedTriEdgeNumber[0])
		k, Kmax, myThread = pm.GetLocalK(int(e.ConnectedTris[0]))
		Np1               = c.dfr.N + 1
		Np12              = float64(Np1 * Np1)
		shift             = edgeNum * Nedge
		edgeLen           = e.GetEdgeLength()
	)
	fs := 0.5 * Np12 * edgeLen / Jdet[myThread].DataP[k]
	edgeMax := -100.
	for i := shift; i < shift+Nedge; i++ {
		ind := k + i*Kmax
		C := c.FS.GetFlowFunction(Q_Face[myThread], ind, SoundSpeed)
		U := c.FS.GetFlowFunction(Q_Face[myThread], ind, Velocity)
		waveSpeed := fs * (U + C)
		waveSpeedMax = math.Max(waveSpeed, waveSpeedMax)
		if waveSpeed > edgeMax {
			edgeMax = waveSpeed
		}
	}
	if edgeMax > DT[myThread].DataP[k] {
		DT[myThread].DataP[k] = edgeMax
	}
	if e.NumConnectedTris == 2 { // Add the wavespeed to the other tri connected to this edge if needed
		k, Kmax, myThread = pm.GetLocalK(int(e.ConnectedTris[1]))
		if edgeMax > DT[myThread].DataP[k] {
			DT[myThread].DataP[k] = edgeMax
		}
	}
	return
}

func (c *Euler) calculateNumericalFlux(Nedge, kL, KmaxL, edgeNumberL, myThreadL, kR, KmaxR, edgeNumberR, myThreadR int,
	normal0 [2]float64, normalFlux, avgSolution [][4]float64,
	Q_Face [][4]utils.Matrix) {
	var (
		shiftL, shiftR = edgeNumberL * Nedge, edgeNumberR * Nedge
	)
	switch c.FluxCalcAlgo {
	case FLUX_Average:
		c.AvgFlux(kL, kR, KmaxL, KmaxR, shiftL, shiftR, Q_Face[myThreadL], Q_Face[myThreadR], normal0, normalFlux)
	case FLUX_LaxFriedrichs:
		c.LaxFlux(kL, kR, KmaxL, KmaxR, shiftL, shiftR, Q_Face[myThreadL], Q_Face[myThreadR], normal0, normalFlux)
	case FLUX_Roe:
		c.RoeFlux(kL, kR, KmaxL, KmaxR, shiftL, shiftR, Q_Face[myThreadL], Q_Face[myThreadR], normal0, normalFlux)
	case FLUX_RoeER:
		c.RoeERFlux(kL, kR, KmaxL, KmaxR, shiftL, shiftR, Q_Face[myThreadL], Q_Face[myThreadR], normal0, normalFlux)
	}
	// Store the average solution for this edge - average of the two connected tri edges
	for i := 0; i < Nedge; i++ {
		indL := kL + (i+shiftL)*KmaxL
		indR := kR + (Nedge-1-i+shiftR)*KmaxR
		for n := 0; n < 4; n++ {
			avgSolution[i][n] = 0.5 * (Q_Face[myThreadL][n].DataP[indL] + Q_Face[myThreadR][n].DataP[indR])
		}
	}
}

//func (c *Euler) calculateBCs(e *DG2D.Edge, k, Kmax, myThread, edgeNumber, Nedge int,
//normal0 [2]float64, Q_Face [][4]utils.Matrix, Time float64, normalFlux, avgSolution [][4]float64) {
func (c *Euler) calculateBCs(e *DG2D.Edge, Nedge int, Time float64,
	k, Kmax, edgeNumber, myThread int,
	normal0 [2]float64, normalFlux, avgSolution [][4]float64,
	Q_Face [][4]utils.Matrix) {
	var (
		calculateNormalFlux bool
		shift               = edgeNumber * Nedge
	)
	calculateNormalFlux = true
	switch e.BCType {
	case types.BC_Far:
		c.FarBC(k, Kmax, shift, Q_Face[myThread], normal0)
	case types.BC_IVortex:
		c.IVortexBC(Time, k, Kmax, shift, Q_Face[myThread], normal0)
	case types.BC_Wall, types.BC_Cyl:
		calculateNormalFlux = false
		c.WallBC(k, Kmax, Q_Face[myThread], shift, normal0, normalFlux) // Calculates normal flux directly
	case types.BC_PeriodicReversed, types.BC_Periodic:
		// One edge of the Periodic BC leads to calculation of both sides within the connected tris section, so noop here
		return
	}
	if calculateNormalFlux {
		var Fx, Fy [4]float64
		for i := 0; i < Nedge; i++ {
			ie := i + shift
			ind := k + ie*Kmax
			Fx, Fy = c.CalculateFlux(Q_Face[myThread], ind)
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
			avgSolution[i][n] = Q_Face[myThread][n].DataP[ind]
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
			ind2 := kGlobal + KmaxGlobal*edgeNum
			IInII := dfr.IInII.DataP[ind2]
			for n := 0; n < 4; n++ {
				nFlux, sign := c.EdgeStore.GetEdgeValues(EdgeFlux, myThread, k, n, edgeNum, dfr)
				rtD := F_RT_DOF[n].DataP
				for i := 0; i < Nedge; i++ {
					// Place normed/scaled flux into the RT element space
					ind := k + (2*Nint+i+shift)*Kmax
					if sign > 0 {
						rtD[ind] = nFlux[i] * IInII
					} else {
						rtD[ind] = -nFlux[Nedge-i-1] * IInII
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

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
	Fluxes       [][4]utils.Matrix
	StorageIndex map[types.EdgeKey]int // Index into normal flux storage using edge key
	PMap         *utils.PartitionMap
	Nedge        int
}

type ValueType uint8

const (
	NumericalFluxForEuler ValueType = iota
	EdgeQValues
	ViscousNormalFlux
)

func (c *Euler) NewEdgeStorage() (nf *EdgeValueStorage) {
	var (
		NumEdges = len(c.DFR.Tris.Edges)
	)
	nf = &EdgeValueStorage{
		StorageIndex: make(map[types.EdgeKey]int),
		PMap:         c.Partitions,
		Nedge:        c.DFR.FluxElement.NpEdge,
		Fluxes:       make([][4]utils.Matrix, int(ViscousNormalFlux)+1),
	}
	// Allocate memory for fluxes
	for i := range nf.Fluxes {
		for n := 0; n < 4; n++ {
			nf.Fluxes[i][n] = utils.NewMatrix(NumEdges, nf.Nedge)
		}
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
		Nedge   = nf.Nedge
		target  = nf.Fluxes[int(valType)][varNum]
	)
	en := dfr.EdgeNumber[kGlobal+Kmax*edgeNum]
	edgeIndex := nf.StorageIndex[en]
	ind := edgeIndex * Nedge
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

func (nf *EdgeValueStorage) PutEdgeValues(en types.EdgeKey, valType ValueType, EdgeValues [][4]float64) {
	var (
		target = nf.Fluxes[int(valType)]
	)
	// Load the normal flux into the global normal flux storage
	edgeIndex := nf.StorageIndex[en]
	for n := 0; n < 4; n++ {
		for i := 0; i < nf.Nedge; i++ {
			ind := i + edgeIndex*nf.Nedge
			target[n].DataP[ind] = EdgeValues[i][n]
		}
	}
}

func (c *Euler) GetFaceNormal(kGlobal, edgeNumber int) (normal [2]float64) {
	var (
		KmaxGlobal = c.DFR.K
		faceInd    = kGlobal + KmaxGlobal*edgeNumber
	)
	normal = [2]float64{c.DFR.FaceNorm[0].DataP[faceInd], c.DFR.FaceNorm[1].DataP[faceInd]}
	return
}

func (c *Euler) StoreDissipationEdgeFlux(epsilon [][]float64,
	edgeKeys EdgeKeySlice, EdgeQ1 [][4]float64) {
	var (
		Nedge       = c.DFR.FluxElement.NpEdge
		viscousFlux = EdgeQ1
		pm          = c.Partitions
		Nint        = c.DFR.FluxElement.NpInt
	)
	for _, en := range edgeKeys {
		e := c.DFR.Tris.Edges[en]
		var (
			ooedgeLength         = 1. / e.GetEdgeLength()
			kLGlobal             = int(e.ConnectedTris[0])
			kL, KmaxL, myThreadL = pm.GetLocalK(int(e.ConnectedTris[0]))
			edgeNumberL          = int(e.ConnectedTriEdgeNumber[0])
			DissXL, DissYL       = c.Dissipation.DissX[myThreadL], c.Dissipation.DissY[myThreadL]
			shiftL               = Nedge * edgeNumberL
			normalL              = c.GetFaceNormal(kLGlobal, edgeNumberL)
			nxL, nyL             = normalL[0], normalL[1]
			epsilonL             = epsilon[myThreadL][kL]

			kR, KmaxR, myThreadR = pm.GetLocalK(int(e.ConnectedTris[1]))
			edgeNumberR          = int(e.ConnectedTriEdgeNumber[1])
			DissXR, DissYR       = c.Dissipation.DissX[myThreadR], c.Dissipation.DissY[myThreadR]
			shiftR               = Nedge * edgeNumberR
			normalR              = c.GetFaceNormal(kLGlobal, edgeNumberL)
			nxR, nyR             = normalR[0], normalR[1]
			epsilonR             = epsilon[myThreadR][kR]
		)
		switch e.NumConnectedTris {
		case 0:
			panic("unable to handle unconnected edges")
		case 1: // Handle edges with only one triangle
			for n := 0; n < 4; n++ {
				for i := 0; i < Nedge; i++ {
					indL := kL + (2*Nint+shiftL+i)*KmaxL
					viscousFlux[i][n] = nxL*DissXL[n].DataP[indL] + nyL*DissYL[n].DataP[indL]
				}
			}
		case 2: // Handle edges with two connected tris - shared faces
			/*
				Numerical viscous flux Fⁿ
						On each interior face, define
						  {λ∇u} = ½ (λ⁻ ∇u⁻ + λ⁺ ∇u⁺),   [u] = u⁻ − u⁺

						and with n the outward normal from the “−” side, the
						SIPG viscous flux is
						  Fⁿ = {λ∇u} · n − σ {λ}/h (u⁻ − u⁺)

						where {λ} = (λ⁻ + λ⁺) / 2.
					σ   = Cp*p*p                  // choose Cp≈1…10, p=poly degree
					h   = faceSize(e)             // characteristic length of the face

						  • The first term uses the normal n once: on the “+”
							side its normal is −n, so it automatically flips sign.
						  • The penalty −σ {λ}/h (u⁻ − u⁺) does not flip sign.

			*/
			for n := 0; n < 4; n++ {
				for i := 0; i < Nedge; i++ {
					indL := kL + (2*Nint+shiftL+i)*KmaxL
					vFL := nxL*DissXL[n].DataP[indL] + nyL*DissYL[n].DataP[indL]
					indR := kR + (2*Nint+shiftR+Nedge-1-i)*KmaxR // Reversed edge - storage is for primary element
					vFR := nxR*DissXR[n].DataP[indR] + nyR*DissYR[n].DataP[indR]
					viscousFlux[i][n] = 0.5 * (vFL + vFR)
				}
			}
			// σ = Cp * p * p // choose Cp≈1…10, p=poly degree
			Cp := 1.
			Omega := Cp * float64(c.DFR.N*c.DFR.N)
			LambdaAvg := 0.5 * (epsilonL + epsilonR)
			for n := 0; n < 4; n++ {
				// Load edge solution values from solution storage
				edgeQL, _ := c.EdgeStore.GetEdgeValues(EdgeQValues, myThreadL,
					kL, n, edgeNumberL, c.DFR)
				edgeQR, _ := c.EdgeStore.GetEdgeValues(EdgeQValues, myThreadR,
					kR, n, edgeNumberR, c.DFR)
				for i := 0; i < Nedge; i++ {
					ii := Nedge - 1 - i
					viscousFlux[i][n] -=
						(Omega * LambdaAvg * ooedgeLength) * (edgeQL[i] - edgeQR[ii])
				}
			}
		}
		// Load the normal flux into the global normal flux storage
		c.EdgeStore.PutEdgeValues(en, ViscousNormalFlux, viscousFlux)
	}
	return
}

func (c *Euler) CalculateEdgeFlux(Time float64, CalculateDT bool, Jdet, DT []utils.Matrix, Q_Face [][4]utils.Matrix,
	Flux_Face [][2][4]utils.Matrix, edgeKeys EdgeKeySlice, EdgeQ1, EdgeQ2 [][4]float64) (waveSpeedMax float64) {
	var (
		Nedge                 = c.DFR.FluxElement.NpEdge
		numericalFluxForEuler = EdgeQ1
		edgeQValues           = EdgeQ2
		pm                    = c.Partitions
	)
	for _, en := range edgeKeys {
		e := c.DFR.Tris.Edges[en]
		var (
			kLGlobal             = int(e.ConnectedTris[0])
			kL, KmaxL, myThreadL = pm.GetLocalK(int(e.ConnectedTris[0]))
			edgeNumberL          = int(e.ConnectedTriEdgeNumber[0])
			normalL              = c.GetFaceNormal(kLGlobal, edgeNumberL)
			shiftL               = edgeNumberL * Nedge
		)
		// Store solution for this edge - use left side only per Cockburn and Shu's algorithm for Laplacian, alternate flux sides
		for i := 0; i < Nedge; i++ {
			indL := kL + (i+shiftL)*KmaxL
			for n := 0; n < 4; n++ {
				edgeQValues[i][n] = Q_Face[myThreadL][n].DataP[indL]
			}
		}
		c.EdgeStore.PutEdgeValues(en, EdgeQValues, edgeQValues)

		for i := range numericalFluxForEuler {
			for n := 0; n < 4; n++ {
				numericalFluxForEuler[i][n] = 0
			}
		}
		switch e.NumConnectedTris {
		case 0:
			panic("unable to handle unconnected edges")
		case 1: // Handle edges with only one triangle - default is edge flux, which will be replaced by a BC flux
			c.calculateNonSharedEdgeFlux(e, Nedge, Time,
				kL, KmaxL, edgeNumberL, myThreadL,
				normalL, numericalFluxForEuler, edgeQValues, Q_Face)
		case 2: // Handle edges with two connected tris - shared faces
			var (
				kR, KmaxR, myThreadR = pm.GetLocalK(int(e.ConnectedTris[1]))
				edgeNumberR          = int(e.ConnectedTriEdgeNumber[1])
			)
			c.calculateSharedEdgeFlux(Nedge, kL, KmaxL, edgeNumberL, myThreadL, kR, KmaxR, edgeNumberR, myThreadR,
				normalL, numericalFluxForEuler, Q_Face, Flux_Face)
		}
		// Load the normal flux into the global normal flux storage
		c.EdgeStore.PutEdgeValues(en, NumericalFluxForEuler, numericalFluxForEuler)

		if CalculateDT {
			waveSpeedMax = c.calculateLocalDT(e, Nedge, Q_Face, Jdet, DT)
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
		Np1               = c.DFR.N + 1
		Np12              = float64(Np1 * Np1)
		shift             = edgeNum * Nedge
		edgeLen           = e.GetEdgeLength()
	)
	fs := 0.5 * Np12 * edgeLen / Jdet[myThread].DataP[k]
	edgeMax := -100.
	for i := shift; i < shift+Nedge; i++ {
		ind := k + i*Kmax
		C := c.FSFar.GetFlowFunction(Q_Face[myThread], ind, SoundSpeed)
		U := c.FSFar.GetFlowFunction(Q_Face[myThread], ind, Velocity)
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

func (c *Euler) calculateSharedEdgeFlux(Nedge, kL, KmaxL, edgeNumberL, myThreadL, kR, KmaxR, edgeNumberR, myThreadR int,
	normalL [2]float64, numericalFluxForEuler [][4]float64, Q_Face [][4]utils.Matrix, Flux_Face [][2][4]utils.Matrix) {
	var (
		shiftL, shiftR = edgeNumberL * Nedge, edgeNumberR * Nedge
		// interpolateFluxNotQ = true
		interpolateFluxNotQ = false
	)
	switch c.FluxCalcAlgo {
	case FLUX_Average:
		c.AvgFlux(kL, kR, KmaxL, KmaxR, shiftL, shiftR, Q_Face[myThreadL], Q_Face[myThreadR], normalL, numericalFluxForEuler)
	case FLUX_LaxFriedrichs:
		if interpolateFluxNotQ {
			c.LaxFlux2(kL, kR, KmaxL, KmaxR, shiftL, shiftR, Q_Face[myThreadL], Q_Face[myThreadR], Flux_Face[myThreadL],
				Flux_Face[myThreadR], normalL, numericalFluxForEuler)
		} else {
			c.LaxFlux(kL, kR, KmaxL, KmaxR, shiftL, shiftR, Q_Face[myThreadL], Q_Face[myThreadR], normalL, numericalFluxForEuler)
		}
	case FLUX_Roe:
		if interpolateFluxNotQ {
			c.RoeFlux2(kL, kR, KmaxL, KmaxR, shiftL, shiftR, Q_Face[myThreadL], Q_Face[myThreadR], Flux_Face[myThreadL],
				Flux_Face[myThreadR], normalL, numericalFluxForEuler)
		} else {
			c.RoeFlux(kL, kR, KmaxL, KmaxR, shiftL, shiftR, Q_Face[myThreadL], Q_Face[myThreadR], normalL, numericalFluxForEuler)
		}
	case FLUX_RoeER:
		c.RoeERFlux(kL, kR, KmaxL, KmaxR, shiftL, shiftR, Q_Face[myThreadL], Q_Face[myThreadR], normalL, numericalFluxForEuler)
	}
}

func (c *Euler) calculateNonSharedEdgeFlux(e *DG2D.Edge, Nedge int, Time float64,
	k, Kmax, edgeNumber, myThread int,
	normal0 [2]float64, numericalFluxForEuler, qFluxForGradient [][4]float64,
	Q_Face [][4]utils.Matrix) {
	var (
		calculateNormalFlux bool
		shift               = edgeNumber * Nedge
	)
	calculateNormalFlux = true
	switch e.BCType {
	case types.BC_None:
		// Do nothing, calculate normal flux
	case types.BC_Far:
		c.FarBC(c.FSFar, k, Kmax, shift, Q_Face[myThread], normal0)
	case types.BC_In:
		c.FarBC(c.FSIn, k, Kmax, shift, Q_Face[myThread], normal0)
	case types.BC_Out:
		c.FarBC(c.FSOut, k, Kmax, shift, Q_Face[myThread], normal0)
	case types.BC_IVortex:
		c.IVortexBC(Time, k, Kmax, shift, Q_Face[myThread], normal0)
	case types.BC_Wall, types.BC_Cyl:
		calculateNormalFlux = false
		c.WallBC(k, Kmax, Q_Face[myThread], shift, normal0, numericalFluxForEuler) // Calculates normal flux directly
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
				numericalFluxForEuler[i][n] = normal0[0]*Fx[n] + normal0[1]*Fy[n]
			}
		}
	}
	return
}

func (c *Euler) SetRTFluxOnEdges(myThread, Kmax int, F_RT_DOF [4]utils.Matrix) {
	var (
		dfr        = c.DFR
		Nedge      = dfr.FluxElement.NpEdge
		Nint       = dfr.FluxElement.NpInt
		KmaxGlobal = c.DFR.K
	)
	for k := 0; k < Kmax; k++ {
		kGlobal := c.Partitions.GetGlobalK(k, myThread)
		for edgeNum := 0; edgeNum < 3; edgeNum++ {
			shift := edgeNum * Nedge
			// nFlux, sign := c.EdgeStore.GetEdgeNormalFlux(kGlobal, edgeNum, DFR)
			ind2 := kGlobal + KmaxGlobal*edgeNum
			IInII := dfr.IInII.DataP[ind2]
			for n := 0; n < 4; n++ {
				nFlux, sign := c.EdgeStore.GetEdgeValues(NumericalFluxForEuler, myThread, k, n, edgeNum, dfr)
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

func (c *Euler) InterpolateSolutionToEdges(Q, Q_Face, Q_Face_P0 [4]utils.Matrix) {
	// Interpolate from solution points to edges using precomputed interpolation matrix
	for n := 0; n < 4; n++ {
		c.DFR.FluxEdgeInterp.Mul(Q[n], Q_Face[n])
		if c.DFR.N > 1 {
			c.DFR.FluxEdgeProject0Interp.Mul(Q[n], Q_Face_P0[n])
		}
	}
	return
}

func (c *Euler) InterpolateSolutionToShockedEdges(sf *DG2D.ModeAliasShockFinder,
	Q_Face, Q_Face_P0 [4]utils.Matrix) {
	var (
		NpEdge      = c.DFR.FluxElement.NpEdge
		NpEdgeTotal = 3 * NpEdge
	)
	for _, k := range sf.ShockCells.Cells() {
		for n := 0; n < 4; n++ {
			for i := 0; i < NpEdgeTotal; i++ {
				Q_Face[n].Set(i, k, Q_Face_P0[n].At(i, k))
			}
		}
	}
	return
}

func (c *Euler) InterpolateSolutionToTargetK(k int,
	Q_Face, Q_Face_P0 [4]utils.Matrix) {
	var (
		NpEdge      = c.DFR.FluxElement.NpEdge
		NpEdgeTotal = 3 * NpEdge
	)
	for n := 0; n < 4; n++ {
		for i := 0; i < NpEdgeTotal; i++ {
			Q_Face[n].Set(i, k, Q_Face_P0[n].At(i, k))
		}
	}
	return
}

func (c *Euler) InterpolateSolutionToTargetEdgeAndK(edge, k int,
	Q_Face, Q_Face_P0 [4]utils.Matrix) {
	var (
		NpEdge = c.DFR.FluxElement.NpEdge
	)
	offset := edge * NpEdge
	for n := 0; n < 4; n++ {
		for i := 0; i < NpEdge; i++ {
			Q_Face[n].Set(i+offset, k, Q_Face_P0[n].At(i+offset, k))
		}
	}
	return
}

func (c *Euler) InterpolateSolutionToEdgesWithEntropyVariables(Q,
	Q_Face [4]utils.Matrix) {
	// Switch to Entropy variables
	SwitchToEntropyVariables(Q, c.FSFar.Gamma)
	// Interpolate from solution points to edges using precomputed interpolation matrix
	for n := 0; n < 4; n++ {
		c.DFR.FluxEdgeInterp.Mul(Q[n], Q_Face[n])
	}
	SwitchToConservedVariables(Q, c.FSFar.Gamma)
	SwitchToConservedVariables(Q_Face, c.FSFar.Gamma)
	return
}

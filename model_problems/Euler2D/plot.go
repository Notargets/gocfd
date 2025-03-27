package Euler2D

import (
	"encoding/binary"
	"os"

	"github.com/notargets/gocfd/DG2D"

	"github.com/notargets/gocfd/types"

	"github.com/notargets/gocfd/utils"
)

func (c *Euler) GetPlotField(Q [4]utils.Matrix,
	plotField FlowFunction) (field utils.Matrix) {
	var (
		Kmax       = c.dfr.K
		Np         = c.dfr.SolutionElement.Np
		NpFlux     = c.dfr.FluxElement.Np
		fld        utils.Matrix
		skipInterp bool
	)
	switch {
	case plotField <= Energy:
		fld = Q[int(plotField)]
	case plotField < ShockFunction:
		fld = utils.NewMatrix(Np, Kmax)
		fldD := fld.DataP
		for ik := 0; ik < Kmax*Np; ik++ {
			fldD[ik] = c.FSFar.GetFlowFunction(Q, ik, plotField)
		}
	case plotField == ShockFunction:
		var sf *DG2D.ModeAliasShockFinder
		if c.Limiter != nil {
			sf = c.Limiter.ShockFinder[0]
		} else {
			sf = c.dfr.NewAliasShockFinder(2)
		}
		fld = utils.NewMatrix(Np, Kmax)
		fldD := fld.DataP
		for k := 0; k < Kmax; k++ {
			qElement := Q[0].Col(k)
			m := sf.ShockIndicator(qElement.DataP)
			for i := 0; i < Np; i++ {
				ik := k + i*Kmax
				fldD[ik] = m
			}
		}
	case plotField == EpsilonDissipation:
		field = c.Dissipation.GetScalarEpsilonPlotField(c)
		skipInterp = true
	case plotField == EpsilonDissipationC0:
		field = c.Dissipation.GetC0EpsilonPlotField(c)
		skipInterp = true
	case plotField == XGradientDensity || plotField == YGradientDensity ||
		plotField == XGradientXMomentum || plotField == YGradientXMomentum ||
		plotField == XGradientYMomentum || plotField == YGradientYMomentum ||
		plotField == XGradientEnergy || plotField == YGradientEnergy:
		GradX, GradY := utils.NewMatrix(NpFlux, Kmax), utils.NewMatrix(NpFlux, Kmax)
		DOFX, DOFY := utils.NewMatrix(NpFlux, Kmax), utils.NewMatrix(NpFlux, Kmax)
		var varNum int
		switch plotField {
		case XGradientDensity, YGradientDensity:
			varNum = 0
		case XGradientXMomentum, YGradientXMomentum:
			varNum = 1
		case XGradientYMomentum, YGradientYMomentum:
			varNum = 2
		case XGradientEnergy, YGradientEnergy:
			varNum = 3
		}
		c.GetSolutionGradientUsingRTElement(-1, varNum, Q, GradX, GradY, DOFX, DOFY)
		if plotField < 300 {
			field = GradX
		} else {
			field = GradY
		}
		skipInterp = true
	}
	if !skipInterp {
		// field = c.dfr.FluxInterp.Mul(fld)
		// fmt.Printf("Function: %s, Min/Max = %.2f/%.2f\n",
		// 	plotField.String(), fld.Min(), fld.Max())
		field = c.dfr.GraphInterp.Mul(fld)
		c.GetFirstOrderEdgeProjection_ForGraphing(fld, &field)
		field = field.Transpose()
	}
	return
}

func (c *Euler) SerializeBCs() (BCXY map[string][][]float32) {
	// We output the XY coordinates of boundary conditions
	var (
		bcEdges = c.dfr.BCEdges
	)
	BCXY = make(map[string][][]float32)
	for _, name := range bcEdges.ListNames() {
		bName := types.BCTAG(name)
		if _, present := bcEdges[bName]; present {
			// We need to (re)construct contiguous BC lists of indices from the edge map
			var vLast int
			ii := -1
			for i, edge := range bcEdges[bName] {
				verts := edge.GetVertices()
				v1, v2 := verts[0], verts[1]
				if i == 0 || v1 != vLast {
					ii++
					BCXY[name] = append(BCXY[name], make([]float32, 1))
					BCXY[name][ii] = append(BCXY[name][ii], float32(c.dfr.VX.DataP[v1]))
					BCXY[name][ii] = append(BCXY[name][ii], float32(c.dfr.VY.DataP[v1]))
					vLast = v1
				}
				BCXY[name][ii] = append(BCXY[name][ii], float32(c.dfr.VX.DataP[v2]))
				BCXY[name][ii] = append(BCXY[name][ii], float32(c.dfr.VY.DataP[v2]))
				vLast = v2
			}
		}
	}
	return
}

func (c *Euler) AppendBCsToMeshFile(fileName string) {
	var (
		err     error
		file    *os.File
		bcEdges = c.dfr.BCEdges
	)
	file, err = os.OpenFile(fileName, 0, os.ModeAppend)
	if err != nil {
		panic(err)
	}
	defer file.Close()

	// We output the XY coordinates of boundary conditions
	var nBCs int64
	for _, name := range bcEdges.ListNames() {
		bName := types.BCTAG(name)
		if _, present := bcEdges[bName]; present {
			nBCs++
		}
	}
	binary.Write(file, binary.LittleEndian, nBCs)
	X := c.dfr.VX
	Y := c.dfr.VY
	var x1, y1, x2, y2 float32
	for _, name := range bcEdges.ListNames() {
		bName := types.BCTAG(name)
		if _, present := bcEdges[bName]; present {
			var fString [16]byte
			copy(fString[:], name)
			binary.Write(file, binary.LittleEndian, fString)
			bcLen := int64(len(bcEdges[bName]))
			binary.Write(file, binary.LittleEndian, bcLen)
			for _, edge := range bcEdges[bName] {
				verts := edge.GetVertices()
				v1, v2 := verts[0], verts[1]
				x1, y1 = float32(X.DataP[v1]), float32(Y.DataP[v1])
				x2, y2 = float32(X.DataP[v2]), float32(Y.DataP[v2])
				binary.Write(file, binary.LittleEndian, x1)
				binary.Write(file, binary.LittleEndian, y1)
				binary.Write(file, binary.LittleEndian, x2)
				binary.Write(file, binary.LittleEndian, y2)
			}
		}
	}
	return
}

func (c *Euler) GetFirstOrderEdgeProjection_ForGraphing(QField utils.Matrix,
	QGraph *utils.Matrix) {
	var (
		dfr         = c.dfr
		NpInt, KMax = QField.Dims()
		NpEdge      = dfr.FluxElement.NpEdge
		NpGraph     = 3*(1+NpEdge) + NpInt
		efi         = dfr.GetEdgeSegmentFluxIndex()
	)
	if QGraph.IsEmpty() {
		QGraphP := utils.NewMatrix(NpGraph, KMax)
		*QGraph = QGraphP
	}
	c.ShockFinder.UpdateShockedCells(QField)
	// d1, d2 := QGraph.Dims()
	for _, k := range c.ShockFinder.ShockCells.Cells() {
		// fmt.Printf("Kmax = %d, k = %d\n", KMax, k)
		// fmt.Printf("QGraph Np/KMax = %d/%d\n", d1, d2)
		// There are NpEdge-1 interior points supporting reconstruction of
		// NpEdge-1 sub-segments on each of the three edges
		// Here we will use the two adjoining corner segments to construct
		// the vertex value and we'll average segments to create interior
		// node values

		// Below in reconstructed efi coordinates:
		// Nseg = NpE-1,  Nefi = 3Nseg
		Nseg := NpEdge - 1
		Nefi := 3 * Nseg
		if len(efi.InteriorPtsIndex) != Nefi {
			panic("edge segment point count is not correct")
		}
		// We'll do the vertices first
		// NpE == NpEdge, below in GraphMesh coordinates
		// v1,   Edge1,   v2,        Edge2,         v3         Edge3
		//  0,  1->NpE,  NpE+1, (NpE+2)->2*NpE+1, 2*NpE+2, 2*NpE+3->3*NpE+2

		// Efi coordinates:
		//      v1              v2                 v3
		// 0+(Nefi-1)/2  ((Nseg-1)+Nseg)/2  (2Nseg-1)+2Nseg/2
		// fmt.Println("NpEdge, Nseg, Nefi, NpInt = ", NpEdge, Nseg, Nefi, NpInt)
		// fmt.Println("InteriorPtsIndex = ", efi.InteriorPtsIndex)
		v1a, v1b := efi.InteriorPtsIndex[0], efi.InteriorPtsIndex[Nefi-1]
		// fmt.Println("v1a/b = ", v1a, v1b)
		QGraph.Set(0, k, 0.5*(QField.At(v1a, k)+QField.At(v1b, k)))
		v2a, v2b := efi.InteriorPtsIndex[Nseg-1], efi.InteriorPtsIndex[Nseg]
		// fmt.Println("v2a/b = ", v2a, v2b)
		QGraph.Set(NpEdge+1, k, 0.5*(QField.At(v2a, k)+QField.At(v2b, k)))
		v3a, v3b := efi.InteriorPtsIndex[2*Nseg-1], efi.InteriorPtsIndex[2*Nseg]
		// fmt.Println("v3a/b = ", v3a, v3b)
		QGraph.Set(2*NpEdge+2, k, 0.5*(QField.At(v3a, k)+QField.At(v3b, k)))

		// Beginning/End Edge Points (0-based), inclusive:
		// Edge1			Edge2				Edge3
		// b:0 e:Nseg-1		b:Nseg e:2*Nseg-1	b:2*Nseg e:3*Nseg-1 or e:Nefi-1

		// Below are loop indices in efi coordinates (non-inclusive)
		// Edge1         Edge2              Edge3
		// 0->Nseg       Nseg->2Nseg        2Nseg->Nefi

		// NpE == NpEdge, below in GraphMesh coordinates (ranges non-inclusive)
		// v1,   Edge1,   v2,        Edge2,         v3         Edge3
		//  0,  1->NpE+1, NpE+1, (NpE+2)->2*NpE+2, 2*NpE+2, 2*NpE+3->3*NpE+3
		var skEdge, skSeg int
		for n := 0; n < 3; n++ { // Each edge
			// Beginning point of range, excluding vertex
			// QGraph.Set(n*NpEdge+n+1, k, QField.At(efi.InteriorPtsIndex[skSeg], k))
			skEdge = 1 + n*(NpEdge+1) // Skip the vertex
			sleft := efi.InteriorPtsIndex[skSeg]
			QGraph.Set(skEdge, k, QField.At(sleft, k))
			skEdge++
			// Interior range
			for i := 1; i < NpEdge-1; i++ { // Averaging segment values
				sright := efi.InteriorPtsIndex[skSeg+1]
				// QGraph.Set(i+n*NpEdge+n+1, k, 0.5*(QField.At(sleft, k)+QField.At(sright,k)))
				QGraph.Set(skEdge, k, 0.5*(QField.At(sleft, k)+QField.At(sright, k)))
				skSeg++
				skEdge++
				sleft = sright
			}
			// End point of range, excluding vertex
			// QGraph.Set((n+1)*NpEdge+n, k, QField.At(efi.InteriorPtsIndex[skSeg], k))
			QGraph.Set(skEdge, k, QField.At(sleft, k))
			skEdge++
			skSeg++
		}
		// fmt.Println("skEdge, skSeg = ", skEdge, skSeg)
	}
}

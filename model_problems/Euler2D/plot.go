package Euler2D

import (
	"encoding/binary"
	"os"

	"github.com/notargets/gocfd/DG2D"

	"github.com/notargets/gocfd/types"

	"github.com/notargets/gocfd/utils"
)

func (c *Euler) GetPlotField(Q [4]utils.Matrix, plotField FlowFunction) (field utils.Matrix) {
	var (
		Kmax       = c.DFR.K
		Np         = c.DFR.SolutionElement.Np
		NpFlux     = c.DFR.FluxElement.Np
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
			sf = c.DFR.NewAliasShockFinder(2)
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
		field = c.DFR.GraphInterp.Mul(fld)
		// c.GetFirstOrderEdgeProjection_ForGraphing(fld, &field)
		field = field.Transpose()
	}
	return
}

func (c *Euler) SerializeBCs() (BCXY map[string][][]float32) {
	// We output the XY coordinates of boundary conditions
	var (
		bcEdges = c.DFR.BCEdges
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
					BCXY[name][ii] = append(BCXY[name][ii], float32(c.DFR.VX.DataP[v1]))
					BCXY[name][ii] = append(BCXY[name][ii], float32(c.DFR.VY.DataP[v1]))
					vLast = v1
				}
				BCXY[name][ii] = append(BCXY[name][ii], float32(c.DFR.VX.DataP[v2]))
				BCXY[name][ii] = append(BCXY[name][ii], float32(c.DFR.VY.DataP[v2]))
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
		bcEdges = c.DFR.BCEdges
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
	X := c.DFR.VX
	Y := c.DFR.VY
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

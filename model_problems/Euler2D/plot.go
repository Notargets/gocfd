package Euler2D

import (
	"encoding/binary"
	"os"

	"github.com/notargets/gocfd/types"

	"github.com/notargets/gocfd/utils"
)

func (c *Euler) GetPlotField(Q [4]utils.Matrix, plotField FlowFunction) (field utils.Matrix) {
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
		var sf *ModeAliasShockFinder
		if c.Limiter != nil {
			sf = c.Limiter.ShockFinder[0]
		} else {
			sf = NewAliasShockFinder(c.dfr.SolutionElement, 2)
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
		field = c.dfr.FluxInterp.Mul(fld)
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

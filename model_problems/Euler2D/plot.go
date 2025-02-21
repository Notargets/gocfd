package Euler2D

import (
	"encoding/binary"
	"fmt"
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

func (c *Euler) SaveOutputMesh(fileName string) {
	var (
		err     error
		file    *os.File
		bcEdges = c.dfr.BCEdges
		gm      = c.dfr.OutputMesh()
	)
	file, err = os.Create(fileName)
	if err != nil {
		panic(err)
	}
	defer file.Close()

	lenTriVerts := int64(3 * len(gm.Triangles))
	triVerts := make([]int64, lenTriVerts)
	for k, tri := range gm.Triangles {
		elp := 3 * k
		triVerts[elp] = int64(tri.Nodes[0])
		triVerts[elp+1] = int64(tri.Nodes[1])
		triVerts[elp+2] = int64(tri.Nodes[2])
	}
	lenXYCoords := int64(len(gm.Geometry))
	xy := make([]float64, lenXYCoords*2) // combined X,Y
	for i, pt := range gm.Geometry {
		xy[2*i] = float64(pt.X[0])
		xy[2*i+1] = float64(pt.X[1])
	}
	fmt.Printf("Number of Coordinate Pairs: %d\n", lenXYCoords)
	fmt.Printf("Number of Original Triangle Elements: %d\n", c.dfr.K)
	fmt.Printf("Number of RT Triangle Elements: %d\n", lenTriVerts/3)
	nDimensions := int64(2) // 2D
	binary.Write(file, binary.LittleEndian, nDimensions)
	binary.Write(file, binary.LittleEndian, lenTriVerts)
	binary.Write(file, binary.LittleEndian, triVerts)
	binary.Write(file, binary.LittleEndian, lenXYCoords)
	binary.Write(file, binary.LittleEndian, xy)
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

func (c *Euler) SavePlotFunction(fI []float32, fileName string, nSteps int) {
	var (
		err error
	)
	if nSteps == 1 {
		c.SolutionOutputFile, err = os.Create(fileName)
		if err != nil {
			panic(err)
		}
		fmt.Printf("Length of scalar output field: %d\n", len(fI))
		fmt.Printf("Solution Element Np = %d\n", c.dfr.SolutionElement.Np)
		fmt.Printf("RT Element Np = %d\n", c.dfr.FluxElement.Np)
	}
	fMin, fMax, fAve := getFRange(fI)
	fmt.Printf("FMin, FMax, FAve: %f, %f, %f\n", fMin, fMax, fAve)
	file := c.SolutionOutputFile
	binary.Write(file, binary.LittleEndian, int64(len(fI)))
	binary.Write(file, binary.LittleEndian, fI)
	return
}

func getFRange(F []float32) (fMin, fMax, fAve float32) {
	var (
		fSum  float32
		count float32
	)
	fMin = F[0]
	fMax = fMin
	fSum = 0
	for _, f := range F {
		if f < fMin {
			fMin = f
		}
		if f > fMax {
			fMax = f
		}
		fSum += f
		count++
	}
	fAve = fSum / count
	return
}

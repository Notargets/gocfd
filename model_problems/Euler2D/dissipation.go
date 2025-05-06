package Euler2D

import (
	"fmt"
	"math"
	"sort"

	"github.com/notargets/gocfd/DG2D"

	"github.com/notargets/gocfd/utils"
)

type VertexToElement [][3]int32 // Vertex id is the first int32, element ID is the next, threadID third

func (ve VertexToElement) Len() int      { return len(ve) }
func (ve VertexToElement) Swap(i, j int) { ve[i], ve[j] = ve[j], ve[i] }
func (ve VertexToElement) Less(i, j int) bool {
	if ve[i][0] < ve[j][0] {
		return true
	} else if ve[i][0] == ve[j][0] {
		if ve[i][1] < ve[j][1] {
			return true
		} else if ve[i][1] == ve[j][1] {
			if ve[i][2] < ve[j][2] {
				return true
			}
		}
	}
	return false
}

func NewVertexToElement(EtoV utils.Matrix) (VtoE VertexToElement) {
	var (
		Kmax, Nverts = EtoV.Dims()
	)
	if Nverts != 3 {
		msg := fmt.Errorf("EtoV should have dimensions [Kmax,3] was [%d,%d]", Kmax, Nverts)
		panic(msg)
	}
	VtoE = make(VertexToElement, Kmax*3)
	var ii int
	for k := 0; k < Kmax; k++ {
		for i := 0; i < 3; i++ {
			VtoE[ii] = [3]int32{int32(EtoV.At(k, i)), int32(k), 0}
			ii++
		}
	}
	sort.Sort(VtoE)
	return
}

func (ve VertexToElement) Shard(pm *utils.PartitionMap) (veSharded []VertexToElement) {
	var (
		NPar             = pm.ParallelDegree
		lve              = len(ve)
		VertexPartitions = utils.NewPartitionMap(NPar, lve) // This has to be re-done to honor vertex grouping
		ib               int
		vNum             int32
	)
	veSharded = make([]VertexToElement, NPar)
	approxBucketSize := VertexPartitions.GetBucketDimension(0)
	getShardedPair := func(vve [3]int32, pm *utils.PartitionMap) (vves [3]int32) {
		nodeIDSharded, _, threadID := pm.GetLocalK(int(vve[1]))
		vves = [3]int32{vve[0], int32(nodeIDSharded), int32(threadID)}
		return
	}
	_ = getShardedPair
	for np := 0; np < NPar; np++ {
		for i := 0; i < approxBucketSize; i++ {
			veSharded[np] = append(veSharded[np], getShardedPair(ve[ib], pm))
			// veSharded[np] = append(veSharded[np], ve[ib])
			ib++
			if ib == lve {
				return
			}
		}
		vNum = ve[ib][0]
		for ib < lve && ve[ib][0] == vNum {
			veSharded[np] = append(veSharded[np], getShardedPair(ve[ib], pm))
			// veSharded[np] = append(veSharded[np], ve[ib])
			ib++
			if ib == lve {
				return
			}
		}
	}
	return
}

type ScalarDissipation struct {
	VtoE                       []VertexToElement   // Sharded vertex to element map, [2] is [vertID, ElementID_Sharded]
	EtoV                       []utils.Matrix      // Sharded Element to Vertex map, Kx3
	Epsilon                    []utils.Matrix      // Sharded Np x Kmax, Interpolated from element vertices
	EpsilonScalar              [][]float64         // Sharded scalar value of dissipation, one per element
	DissDOF, DissDOF2, DissDiv []utils.Matrix      // Sharded NpFlux x Kmax, DOF for Gradient calculation using RT
	DissX, DissY               [][4]utils.Matrix   // Sharded NpFlux x Kmax, R and S derivative of dissipation field
	EpsVertex                  []float64           // NVerts x 1, Aggregated (Max) of epsilon surrounding each vertex, Not sharded
	PMap                       *utils.PartitionMap // Partition map for the solution shards in K
	Clipper                    utils.Matrix        // Matrix used to clip the topmost mode from the solution polynomial, used in shockfinder
	dfr                        *DG2D.DFR2D
	ShockFinder                []*DG2D.ModeAliasShockFinder
	Kappa                      float64
	BaryCentricCoords          utils.Matrix // A thruple(lam0,lam1,lam2) for interpolation for each interior point, Npx3
	VertexEpsilonValues        []utils.Matrix
}

func NewScalarDissipation(kappa float64, dfr *DG2D.DFR2D, pm *utils.PartitionMap) (sd *ScalarDissipation) {
	var (
		NPar   = pm.ParallelDegree
		el     = dfr.SolutionElement
		Np     = el.Np
		NpFlux = dfr.FluxElement.Np
		order  = float64(el.N)
		NVerts = dfr.VX.Len()
	)
	_ = order
	sd = &ScalarDissipation{
		EpsVertex:           make([]float64, NVerts),
		EpsilonScalar:       make([][]float64, NPar),    // Viscosity, constant over the element
		Epsilon:             make([]utils.Matrix, NPar), // Epsilon field, expressed over solution points
		VertexEpsilonValues: make([]utils.Matrix, NPar), // Epsilon field, expressed over solution points
		DissDOF:             make([]utils.Matrix, NPar),
		DissDOF2:            make([]utils.Matrix, NPar),
		DissDiv:             make([]utils.Matrix, NPar),
		DissX:               make([][4]utils.Matrix, NPar),
		DissY:               make([][4]utils.Matrix, NPar),
		VtoE:                NewVertexToElement(dfr.Tris.EToV).Shard(pm),
		PMap:                pm,
		dfr:                 dfr,
		ShockFinder:         make([]*DG2D.ModeAliasShockFinder, NPar),
		Kappa:               5.,
	}
	sd.EtoV = sd.shardEtoV(dfr.Tris.EToV)
	sd.createInterpolationStencil()
	if kappa != 0. {
		sd.Kappa = kappa
	}
	for np := 0; np < NPar; np++ {
		Kmax := pm.GetBucketDimension(np)
		sd.Epsilon[np] = utils.NewMatrix(NpFlux, Kmax)
		sd.VertexEpsilonValues[np] = utils.NewMatrix(3, Kmax)
		sd.EpsilonScalar[np] = make([]float64, Kmax)
		sd.DissDOF[np] = utils.NewMatrix(NpFlux, Kmax)
		sd.DissDOF2[np] = utils.NewMatrix(NpFlux, Kmax)
		sd.DissDiv[np] = utils.NewMatrix(Np, Kmax)
		for n := 0; n < 4; n++ {
			sd.DissX[np][n] = utils.NewMatrix(NpFlux, Kmax)
			sd.DissY[np][n] = utils.NewMatrix(NpFlux, Kmax)
		}
		sd.ShockFinder[np] = dfr.NewAliasShockFinder(sd.Kappa)
	}
	/*
		The "Clipper" matrix drops the last mode from the polynomial and forms an alternative field of values at the node
		points based on a polynomial with one less term. In other words, if we have a polynomial of degree "p", expressed
		as values at Np node points, multiplying the Node point values vector by Clipper produces an alternative version
		of the node values based on truncating the last polynomial mode.
	*/
	sd.Clipper = el.JB2D.V.Mul(dfr.CutoffFilter2D(el.N, el.N-1,
		0)).Mul(el.JB2D.Vinv)
	return
}

func (sd *ScalarDissipation) shardEtoV(EtoV utils.Matrix) (ev []utils.Matrix) {
	var (
		pm = sd.PMap
		NP = pm.ParallelDegree
		// KMax, _ = EtoV.Dims()
	)
	ev = make([]utils.Matrix, NP)
	for np := 0; np < NP; np++ {
		KMax := pm.GetBucketDimension(np)
		ev[np] = utils.NewMatrix(KMax, 3)
		kmin, kmax := pm.GetBucketRange(np)
		var klocal int
		for k := kmin; k < kmax; k++ {
			ev[np].Set(klocal, 0, EtoV.At(k, 0))
			ev[np].Set(klocal, 1, EtoV.At(k, 1))
			ev[np].Set(klocal, 2, EtoV.At(k, 2))
			klocal++
		}
		if klocal != KMax {
			msg := fmt.Errorf("dimension incorrect, should be %d, is %d", KMax, klocal)
			panic(msg)
		}
	}
	return
}

func (sd *ScalarDissipation) propagateEpsilonMaxToVertices(myThread int) {
	var (
		VtoE = sd.VtoE[myThread]
	)
	oldVert := -1
	for _, val := range VtoE {
		vert, k, threadID := int(val[0]), int(val[1]), int(val[2])
		if oldVert == vert { // we're in the middle of processing this vert, update normally
			sd.EpsVertex[vert] = max(sd.EpsVertex[vert], sd.EpsilonScalar[threadID][k])
		} else { // we're on a new vertex, reset the vertex value
			sd.EpsVertex[vert] = sd.EpsilonScalar[threadID][k]
			oldVert = vert
		}
	}
}

type ContinuityLevel uint8

const (
	No ContinuityLevel = iota
	C0
	C1
)

func (sd *ScalarDissipation) CalculateEpsilonGradient(c *Euler, cont ContinuityLevel, myThread int, Q [4]utils.Matrix) {
	var (
		dfr    = sd.dfr
		Kmax   = sd.PMap.GetBucketDimension(myThread)
		NpFlux = dfr.FluxElement.Np

		EpsilonScalar       = sd.EpsilonScalar[myThread]
		Epsilon             = sd.Epsilon[myThread]
		DOF, DOF2           = sd.DissDOF[myThread], sd.DissDOF2[myThread]
		DissX, DissY        = sd.DissX[myThread], sd.DissY[myThread]
		EtoV                = sd.EtoV[myThread]
		VertexEpsilonValues = sd.VertexEpsilonValues[myThread]
	)
	if cont == C0 {
		// Interpolate epsilon within each element
		for k := 0; k < Kmax; k++ {
			tri := EtoV.DataP[3*k : 3*k+3]
			v := [3]int{int(tri[0]), int(tri[1]), int(tri[2])}
			for vert := 0; vert < 3; vert++ {
				ind := k + vert*Kmax
				VertexEpsilonValues.DataP[ind] = sd.EpsVertex[v[vert]]
			}
		}
		sd.BaryCentricCoords.Mul(VertexEpsilonValues, Epsilon)
	}
	for n := 0; n < 4; n++ {
		c.GetSolutionGradientUsingRTElement(myThread, n, Q, DissX[n], DissY[n], DOF, DOF2)
		switch cont {
		case No:
			for k := 0; k < Kmax; k++ {
				for i := 0; i < NpFlux; i++ {
					ind := k + Kmax*i
					DissX[n].DataP[ind] *= EpsilonScalar[k] // Scalar viscosity, constant within each k'th element
					DissY[n].DataP[ind] *= EpsilonScalar[k]
				}
			}
		case C0:
			DissX[n].ElMul(Epsilon)
			DissY[n].ElMul(Epsilon)
		}
	}
}

func (sd *ScalarDissipation) AddDissipation(c *Euler, myThread int, Jinv, Jdet utils.Matrix, RHSQ [4]utils.Matrix) {
	/*
		The dissipation term is in the form:
		diss = epsilon*Grad(U)

		dU/dT = -Div(Flux) + Div(diss)
		dU/dT = -Div(Flux) + Div(epsilon*Grad(U))
		dU/dT = -(Div(Flux) - Div(epsilon*Grad(U)))
		dU/dT = -Div(Flux - epsilon*Grad(U))
	*/
	var (
		dfr          = sd.dfr
		Kmax         = sd.PMap.GetBucketDimension(myThread)
		NpInt        = dfr.FluxElement.NpInt
		KmaxGlobal   = sd.PMap.MaxIndex
		DOF          = sd.DissDOF[myThread]
		DIV          = sd.DissDiv[myThread]
		DissX, DissY = sd.DissX[myThread], sd.DissY[myThread]
	)
	for n := 0; n < 4; n++ {
		/*
			Add the DissX and DissY to the RT_DOF using the contravariant transform for the interior
			and IInII for the edges
		*/
		var (
			DiXd, DiYd = DissX[n].DataP, DissY[n].DataP
			NpEdge     = dfr.FluxElement.NpEdge
			DOFd       = DOF.DataP
		)
		for k := 0; k < Kmax; k++ {
			var (
				JdetD   = Jdet.DataP[k]
				JinvD   = Jinv.DataP[4*k : 4*(k+1)]
				IInIId  = dfr.IInII.DataP
				kGlobal = sd.PMap.GetGlobalK(k, myThread)
			)
			for i := 0; i < NpInt; i++ {
				ind := k + Kmax*i
				ind2 := k + Kmax*(i+NpInt)
				DOFd[ind] = JdetD * (JinvD[0]*DiXd[ind] + JinvD[1]*DiYd[ind])
				DOFd[ind2] = JdetD * (JinvD[2]*DiXd[ind] + JinvD[3]*DiYd[ind])
			}
			for edgeNum := 0; edgeNum < 3; edgeNum++ {
				var (
					fInd  = kGlobal + KmaxGlobal*edgeNum
					IInII = IInIId[fInd]
					shift = NpEdge * edgeNum
				)
				edgeFlux, sign := c.EdgeStore.GetEdgeValues(ViscousNormalFlux, myThread, k, n, edgeNum, c.DFR)
				var ii int
				for i := 0; i < NpEdge; i++ {
					ind := k + (2*NpInt+i+shift)*Kmax
					if sign == -1 {
						ii = NpEdge - 1 - i
					} else {
						ii = i
					}
					DOFd[ind] = edgeFlux[ii] * IInII * float64(sign)
				}
			}
		}
		sd.dfr.FluxElement.DivInt.Mul(DOF, DIV)
		for k := 0; k < Kmax; k++ {
			var (
				oojd = 1. / Jdet.DataP[k]
			)
			for i := 0; i < NpInt; i++ {
				ind := k + i*Kmax
				RHSQ[n].DataP[ind] += oojd * DIV.DataP[ind]
			}
		}
	}
}

func (sd *ScalarDissipation) linearInterpolateEpsilon(myThread int) {
	var (
		dfr          = sd.dfr
		Epsilon      = sd.Epsilon[myThread]
		EtoV         = sd.EtoV[myThread]
		NpFlux, KMax = dfr.FluxElement.Np, sd.PMap.GetBucketDimension(myThread)
		R, S         = dfr.FluxElement.R, dfr.FluxElement.S
	)
	vertLinear := func(r, s float64, f [3]float64) (fi float64) {
		var (
			rLen, sLen     = 2., 2.
			drFrac, dsFrac = (r - (-1.)) / rLen, (s - (-1.)) / sLen
			dr, ds         = f[1] - f[0], f[2] - f[0]
		)
		fi = dr*drFrac + ds*dsFrac + f[0]
		return
	}
	// Interpolate epsilon within each element
	for k := 0; k < KMax; k++ {
		tri := EtoV.Row(k).DataP
		v := [3]int{int(tri[0]), int(tri[1]), int(tri[2])}
		eps := [3]float64{sd.EpsVertex[v[0]], sd.EpsVertex[v[1]], sd.EpsVertex[v[2]]}
		for i := 0; i < NpFlux; i++ {
			ind := k + KMax*i
			Epsilon.DataP[ind] = vertLinear(R.DataP[i], S.DataP[i], eps)
		}
	}
}

func (sd *ScalarDissipation) baryCentricInterpolateEpsilon(myThread int) {
	var (
		dfr      = sd.dfr
		Np, KMax = dfr.SolutionElement.Np, sd.PMap.GetBucketDimension(myThread)
		Epsilon  = sd.Epsilon[myThread]
		EtoV     = sd.EtoV[myThread]
	)
	// Interpolate epsilon within each element
	for k := 0; k < KMax; k++ {
		tri := EtoV.Row(k).DataP
		v := [3]int{int(tri[0]), int(tri[1]), int(tri[2])}
		eps := [3]float64{sd.EpsVertex[v[0]], sd.EpsVertex[v[1]], sd.EpsVertex[v[2]]}
		for i := 0; i < Np; i++ {
			ind := k + KMax*i
			bcc := sd.BaryCentricCoords.Row(i).DataP
			Epsilon.DataP[ind] = bcc[0]*eps[0] + bcc[1]*eps[1] + bcc[2]*eps[2]
		}
	}
}

func (sd *ScalarDissipation) GetScalarEpsilonPlotField(c *Euler) (fld utils.Matrix) {
	for np := 0; np < sd.PMap.ParallelDegree; np++ {
		Np, KMax := sd.Epsilon[np].Dims()
		for k := 0; k < KMax; k++ {
			epsK := sd.EpsilonScalar[np][k]
			for i := 0; i < Np; i++ {
				ind := k + KMax*i
				sd.Epsilon[np].DataP[ind] = epsK
			}
		}
	}
	c.RecombineShardsK(sd.Epsilon, &fld)
	return
}

func (sd *ScalarDissipation) GetC0EpsilonPlotField(c *Euler) (fld utils.Matrix) {
	c.RecombineShardsK(sd.Epsilon, &fld)
	return
}

func (sd *ScalarDissipation) CalculateElementViscosity(myThread int, Qall [][4]utils.Matrix) {
	var (
		dfr  = sd.dfr
		Rho  = Qall[myThread][0]
		Eps  = sd.EpsilonScalar[myThread]
		Kmax = sd.PMap.GetBucketDimension(myThread)
		// U          = sd.U[myThread]
		// UClipped   = sd.UClipped[myThread]
		KMaxGlobal = sd.PMap.MaxIndex
		Order      = float64(sd.dfr.N)
		sf         = sd.ShockFinder[myThread]
	)
	/*
		Eps0 wants to be (h/p) and is supposed to be proportional to cell width
		Something like this for the "h" quantity seems right
			Np1  = c.DFR.N + 1
			Np12 = float64(Np1 * Np1)
			edgeLen     = e.GetEdgeLength()
			fs := 0.5 * Np12 * edgeLen / Jdet[bn].DataP[k]
	*/
	for k := 0; k < Kmax; k++ {
		// Get edges for this element
		kGlobal := sd.PMap.GetGlobalK(k, myThread)
		var maxEdgeLen float64
		maxEdgeLen = -1
		for edgeNum := 0; edgeNum < 3; edgeNum++ {
			ind := kGlobal + KMaxGlobal*edgeNum
			edgeLen := dfr.IInII.DataP[ind]
			if edgeLen > maxEdgeLen {
				maxEdgeLen = edgeLen
			}
		}
		eps0 := 0.75 * maxEdgeLen / Order

		for i := 0; i < sd.dfr.SolutionElement.Np; i++ {
			ind := k + i*Kmax
			sf.Qalt.DataP[i] = Rho.DataP[ind]
		}
		sigma := sf.ShockIndicator(sf.Qalt.DataP)
		Eps[k] = eps0 * sigma
	}
}

func (sd *ScalarDissipation) createInterpolationStencil() {
	var (
		Np    = sd.dfr.FluxElement.Np
		R, S  = sd.dfr.FluxElement.R, sd.dfr.FluxElement.S
		RRinv utils.Matrix
		err   error
	)
	sd.BaryCentricCoords = utils.NewMatrix(Np, 3)
	// Set up unit triangle matrix with vertices in order
	RR := utils.NewMatrix(3, 3)
	RR.Set(0, 0, 1.)
	RR.Set(0, 1, 1.)
	RR.Set(0, 2, 1.)

	RR.Set(1, 0, -1)
	RR.Set(2, 0, -1)

	RR.Set(1, 1, 1)
	RR.Set(2, 1, -1)

	RR.Set(1, 2, -1)
	RR.Set(2, 2, 1)
	if RRinv, err = RR.Inverse(); err != nil {
		panic(err)
	}
	C := utils.NewMatrix(3, 1)
	C.DataP[0] = 1
	for i := 0; i < Np; i++ {
		C.DataP[1] = R.DataP[i]
		C.DataP[2] = S.DataP[i]
		LAM := RRinv.Mul(C)
		for ii := 0; ii < 3; ii++ {
			sd.BaryCentricCoords.DataP[ii+3*i] = LAM.DataP[ii]
		}
	}
}

func (c *Euler) modulateQInterp(Q [4]utils.Matrix, QInterp [4]utils.Matrix, sf *DG2D.ModeAliasShockFinder) {
	var (
		// Beta = 5e-2
		Beta = 2.
		// Beta = 6.6
		// Beta        = 10.
		NpInterp, _ = QInterp[0].Dims()
	)
	sf.UpdateShockedCells(Q[0])
	for _, k := range sf.ShockCells.Cells() {
		sigma := c.ShockFinder.GetShockIndicator(Q[0], k)
		alpha := math.Exp(-Beta * sigma)
		// fmt.Printf("sigma, alpha[%d] = %.1f, %.1f\n", k, sigma, alpha)
		uMean := c.getElementMean(Q, k)
		for n := 0; n < 4; n++ {
			for i := 0; i < NpInterp; i++ {
				QInterp[n].Set(i, k, alpha*QInterp[n].At(i, k)+(1.-alpha)*uMean[n])
			}
		}
	}
	return
}

func (c *Euler) getElementMean(Q [4]utils.Matrix, k int) (Umean [4]float64) {
	var (
		N     = c.DFR.SolutionElement.N
		Np, _ = Q[0].Dims()
	)
	bcn := DG2D.WilliamsShunnJamesonCoords[N]
	for n := 0; n < 4; n++ {
		Umean[n] = 0.0
		for i := 0; i < Np; i++ {
			Umean[n] += Q[n].At(i, k) * bcn[i].W
		}
	}
	return
}

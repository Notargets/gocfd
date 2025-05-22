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
	Continuity                 ContinuityLevel
	VtoE                       []VertexToElement   // Sharded vertex to element map, [2] is [vertID, ElementID_Sharded]
	EtoV                       []utils.Matrix      // Sharded Element to Vertex map, Kx3
	EpsilonScalar              []utils.Vector      // Sharded scalar value of dissipation, one per element
	EpsVertex                  []float64           // NVerts x 1, Aggregated (Max) of epsilon surrounding each vertex, Not sharded
	Epsilon                    []utils.Matrix      // Sharded Np x Kmax, Interpolated from element vertices
	SigmaScalar                []utils.Vector      // Sharded scalar value of Shock Sensor, one per element
	SigmaVertex                []float64           // NVerts x 1, // Aggregated (Max) of Sigma surrounding each vertex, Not sharded
	Sigma                      []utils.Matrix      // Sharded Np x Kmax, Interpolated from element vertices
	DissDOF, DissDOF2, DissDiv []utils.Matrix      // Sharded NpFlux x Kmax, DOF for Gradient calculation using RT
	DissX, DissY               [][4]utils.Matrix   // Sharded NpFlux x Kmax, R and S derivative of dissipation field
	PMap                       *utils.PartitionMap // Partition map for the solution shards in K
	dfr                        *DG2D.DFR2D
	ShockFinder                []*DG2D.ModeAliasShockFinder
	Kappa                      float64
	BaryCentricCoords          utils.Matrix // A thruple(lam0,lam1,lam2) for interpolation for each interior point, Npx3
	VertexScratchSpace         []utils.Matrix
	Eps0                       float64 // Constant multiplier for the viscosity
}

func NewScalarDissipation(kappa float64, cont ContinuityLevel, dfr *DG2D.DFR2D,
	pm *utils.PartitionMap) (sd *ScalarDissipation) {
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
		Continuity:         cont,
		EpsilonScalar:      make([]utils.Vector, NPar), // Viscosity, const over the element
		EpsVertex:          make([]float64, NVerts),    // Vertex values of epsilon
		Epsilon:            make([]utils.Matrix, NPar), // Epsilon field, expressed over solution points
		SigmaScalar:        make([]utils.Vector, NPar), // Shock Finder Result, const over element
		SigmaVertex:        make([]float64, NVerts),    // Vertex values of sigma
		Sigma:              make([]utils.Matrix, NPar), // Sigma field, expressed over solution points
		VertexScratchSpace: make([]utils.Matrix, NPar), // Used to interpolate Vertex values to interior
		DissDOF:            make([]utils.Matrix, NPar),
		DissDOF2:           make([]utils.Matrix, NPar),
		DissDiv:            make([]utils.Matrix, NPar),
		DissX:              make([][4]utils.Matrix, NPar),
		DissY:              make([][4]utils.Matrix, NPar),
		VtoE:               NewVertexToElement(dfr.Tris.EToV).Shard(pm),
		PMap:               pm,
		dfr:                dfr,
		ShockFinder:        make([]*DG2D.ModeAliasShockFinder, NPar),
		Kappa:              5.,
	}
	sd.Eps0 = sd.Kappa / 1.5
	sd.EtoV = sd.shardEtoV(dfr.Tris.EToV)
	sd.createInterpolationStencil()
	if kappa != 0. {
		sd.Kappa = kappa
	}
	for np := 0; np < NPar; np++ {
		KMax := pm.GetBucketDimension(np)
		sd.EpsilonScalar[np] = utils.NewVector(KMax)
		sd.Epsilon[np] = utils.NewMatrix(NpFlux, KMax)
		sd.SigmaScalar[np] = utils.NewVector(KMax)
		sd.Sigma[np] = utils.NewMatrix(NpFlux, KMax)
		sd.VertexScratchSpace[np] = utils.NewMatrix(3, KMax)
		sd.DissDOF[np] = utils.NewMatrix(NpFlux, KMax)
		sd.DissDOF2[np] = utils.NewMatrix(NpFlux, KMax)
		sd.DissDiv[np] = utils.NewMatrix(Np, KMax)
		for n := 0; n < 4; n++ {
			sd.DissX[np][n] = utils.NewMatrix(NpFlux, KMax)
			sd.DissY[np][n] = utils.NewMatrix(NpFlux, KMax)
		}
		sd.ShockFinder[np] = dfr.NewAliasShockFinder(sd.Kappa)
	}
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

func (sd *ScalarDissipation) EpsilonSigmaMaxToVertices(myThread int) {
	var (
		VtoE = sd.VtoE[myThread]
	)
	oldVert := -1
	for _, val := range VtoE {
		vert, k, threadID := int(val[0]), int(val[1]), int(val[2])
		if oldVert == vert { // we're in the middle of processing this vert, update normally
			sd.EpsVertex[vert] = max(sd.EpsVertex[vert], sd.EpsilonScalar[threadID].AtVec(k))
			sd.SigmaVertex[vert] = max(sd.SigmaVertex[vert], sd.SigmaScalar[threadID].AtVec(k))
		} else { // we're on a new vertex, reset the vertex value
			sd.EpsVertex[vert] = sd.EpsilonScalar[threadID].AtVec(k)
			sd.SigmaVertex[vert] = sd.SigmaScalar[threadID].AtVec(k)
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

func (sd *ScalarDissipation) InterpolateEpsilonSigma(myThread int,
	EtoV utils.Matrix) {
	var (
		Kmax = sd.PMap.GetBucketDimension(myThread)
	)
	interpolate := func(VertexValuesM, Interior utils.Matrix, VertexValues []float64) {
		// Interpolate epsilon within each element
		for k := 0; k < Kmax; k++ {
			tri := EtoV.DataP[3*k : 3*k+3]
			v := [3]int{int(tri[0]), int(tri[1]), int(tri[2])}
			for vert := 0; vert < 3; vert++ {
				ind := k + vert*Kmax
				VertexValuesM.DataP[ind] = VertexValues[v[vert]]
			}
		}
		sd.BaryCentricCoords.Mul(VertexValuesM, Interior)
		return
	}
	interpolate(sd.VertexScratchSpace[myThread], sd.Epsilon[myThread],
		sd.EpsVertex)
	interpolate(sd.VertexScratchSpace[myThread], sd.Sigma[myThread],
		sd.SigmaVertex)
	return
}

func (sd *ScalarDissipation) CalculateEpsilonGradient(c *Euler, myThread int, Q [4]utils.Matrix) {
	var (
		dfr    = sd.dfr
		Kmax   = sd.PMap.GetBucketDimension(myThread)
		NpFlux = dfr.FluxElement.Np

		EpsilonScalar = sd.EpsilonScalar[myThread]
		Epsilon       = sd.Epsilon[myThread]
		DOF, DOF2     = sd.DissDOF[myThread], sd.DissDOF2[myThread]
		DissX, DissY  = sd.DissX[myThread], sd.DissY[myThread]
	)
	for n := 0; n < 4; n++ {
		c.GetSolutionGradientUsingRTElement(myThread, n, Q, DissX[n], DissY[n], DOF, DOF2)
		switch sd.Continuity {
		case No:
			for k := 0; k < Kmax; k++ {
				epsScalar := EpsilonScalar.AtVec(k)
				for i := 0; i < NpFlux; i++ {
					ind := k + Kmax*i
					DissX[n].DataP[ind] *= epsScalar
					DissY[n].DataP[ind] *= epsScalar
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

func (sd *ScalarDissipation) GetScalarEpsilonPlotField(c *Euler) (fld utils.Matrix) {
	for np := 0; np < sd.PMap.ParallelDegree; np++ {
		Np, KMax := sd.Epsilon[np].Dims()
		for k := 0; k < KMax; k++ {
			epsK := sd.EpsilonScalar[np].AtVec(k)
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

func (sd *ScalarDissipation) CalculateElementViscosity(myThread int) {
	var (
		// dfr         = sd.dfr
		KMax        = sd.PMap.GetBucketDimension(myThread)
		Eps         = sd.EpsilonScalar[myThread]
		SigmaScalar = sd.SigmaScalar[myThread]
	)
	/*
		Eps0 wants to be (h/p) and is supposed to be proportional to cell width
	*/
	for k := 0; k < KMax; k++ {
		kGlobal := sd.PMap.GetGlobalK(k, myThread)
		// hK := GetHK(dfr.EdgeLenMinR.AtVec(kGlobal))
		Eps.Set(k, sd.Eps0*sd.dfr.GetHkVisc(kGlobal)*SigmaScalar.AtVec(k))
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
	// SetScalar up unit triangle matrix with vertices in order
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
		NpInterp, KMax = QInterp[0].Dims()
	)
	sf.UpdateShockedCells(Q[0])
	for _, k := range sf.ShockCells.Cells() {
		sigma := c.ShockFinder.GetShockIndicator(Q[0], k)
		alpha := math.Exp(-Beta * sigma)
		// fmt.Printf("sigma, alpha[%d] = %.1f, %.1f\n", k, sigma, alpha)
		uMean := c.getElementMean(Q, k)
		for n := 0; n < 4; n++ {
			for i := 0; i < NpInterp; i++ {
				ind := k + i*KMax
				QInterp[n].DataP[ind] = alpha*QInterp[n].DataP[ind] +
					(1.-alpha)*uMean[n]
			}
		}
	}
	return
}

func (c *Euler) getElementMean(Q [4]utils.Matrix, k int) (Umean [4]float64) {
	var (
		N        = c.DFR.SolutionElement.N
		Np, KMax = Q[0].Dims()
	)
	bcn := DG2D.WilliamsShunnJamesonCoords[N]
	for n := 0; n < 4; n++ {
		Umean[n] = 0.0
		for i := 0; i < Np; i++ {
			Umean[n] += Q[n].DataP[k+i*KMax] * bcn[i].W
		}
	}
	return
}

func (sd *ScalarDissipation) UpdateShockFinderSigma(myThread int, Se utils.Vector) {
	/*
		Original method by Persson, constants chosen to match Zhiqiang, et. al.
	*/
	var (
		SigmaScalar = sd.SigmaScalar[myThread]
		kappa       = sd.Kappa
		N           = sd.dfr.SolutionElement.N
		KMax, _     = Se.Dims()
		S0          = 4. / math.Pow(float64(N+1), 4.)
		left, right = S0 - kappa, S0 + kappa
		ookappa     = 0.5 / kappa
		sigma       float64
	)
	for k := 0; k < KMax; k++ {
		se := Se.AtVec(k)
		switch {
		case se < left:
			sigma = 0.
		case se >= left && se <= right:
			sigma = 0.5 * (1. + math.Sin(math.Pi*ookappa*(se-S0)))
		case se > right:
			sigma = 1.
		}
		SigmaScalar.Set(k, sigma)
	}
	return
}

func (sd *ScalarDissipation) LimitSolution(myThread int, Q [4]utils.Matrix,
	QMean [4]utils.Vector) {
	// Note that this approach is equivalent to applying the limiter to modes
	// 1 and higher of the polynomial for the element,
	// as the mean is actually the mode1 value for the polynomial when we
	// have an orthogonal basis. By computing the mean and doing it this way,
	// it's faster to compute.
	var (
		Np, Kmax    = Q[0].Dims()
		SigmaScalar = sd.SigmaScalar[myThread]
		// Sigma       = sd.Sigma[myThread]
	)
	for n := 0; n < 4; n++ {
		for k := 0; k < Kmax; k++ {
			// Smooth ramp accelerator 0-1 for sigma
			alpha := math.Sin(0.5 * math.Pi * SigmaScalar.AtVec(k))
			qmean := QMean[n].AtVec(k)
			for i := 0; i < Np; i++ {
				ind := k + Kmax*i
				q := Q[n].DataP[ind]
				Q[n].DataP[ind] = (1.-alpha)*q + alpha*qmean
			}
		}
	}
	return
}

func (sd *ScalarDissipation) LimitFilterSolution(myThread int,
	Q [4]utils.Matrix, QScratch utils.Matrix, ShockSensor *DG2D.ModeAliasShockFinder) {
	// Note that this approach is equivalent to applying the limiter to modes
	// 1 and higher of the polynomial for the element,
	// as the mean is actually the mode1 value for the polynomial when we
	// have an orthogonal basis. By computing the mean and doing it this way,
	// it's faster to compute.
	var (
		Np, Kmax    = Q[0].Dims()
		SigmaScalar = sd.SigmaScalar[myThread]
		mf          = ShockSensor.ModeFilter
		el          = sd.dfr.SolutionElement
		jb2d        = el.JB2D
		V           = jb2d.V
		Vinv        = jb2d.Vinv
	)
	_ = mf
	for n := 0; n < 4; n++ {
		// Transform Q into modal
		Vinv.Mul(Q[n], QScratch)
		for i := 1; i < Np; i++ {
			for k := 0; k < Kmax; k++ {
				sigma := SigmaScalar.AtVec(k)
				alpha := math.Sin(0.5 * math.Pi * sigma)
				// alpha := math.Sqrt(math.Sin(0.5 * math.Pi * sigma))
				ind := k + Kmax*i
				// QScratch.DataP[ind] *= mf[i] * (1. - alpha)
				QScratch.DataP[ind] *= 1. - alpha
			}
		}
		// Transform QScratch into Q nodal
		V.Mul(QScratch, Q[n])
	}
	return
}

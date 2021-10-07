package Euler2D

import (
	"fmt"
	"math"
	"sort"
	"sync"

	"github.com/notargets/gocfd/DG2D"

	"github.com/notargets/gocfd/utils"
)

type VertexToElement [][3]int32 // Vertex id is the first int32, element ID is the next, threadID third

func (ve VertexToElement) Len() int           { return len(ve) }
func (ve VertexToElement) Swap(i, j int)      { ve[i], ve[j] = ve[j], ve[i] }
func (ve VertexToElement) Less(i, j int) bool { return ve[i][0] < ve[j][0] }
func (ve VertexToElement) Sort()              { sort.Sort(ve) }

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
	VtoE.Sort()
	return
}

func (ve VertexToElement) Shard(pm *PartitionMap) (veSharded []VertexToElement) {
	var (
		NPar             = pm.ParallelDegree
		lve              = len(ve)
		VertexPartitions = NewPartitionMap(NPar, lve) // This has to be re-done to honor vertex grouping
		ib               int
		vNum             int32
	)
	veSharded = make([]VertexToElement, NPar)
	approxBucketSize := VertexPartitions.GetBucketDimension(0)
	getShardedPair := func(vve [3]int32, pm *PartitionMap) (vves [3]int32) {
		nodeIDSharded, _, threadID := pm.GetLocalK(int(vve[1]))
		vves = [3]int32{vve[0], int32(nodeIDSharded), int32(threadID)}
		return
	}
	_ = getShardedPair
	for np := 0; np < NPar; np++ {
		for i := 0; i < approxBucketSize; i++ {
			veSharded[np] = append(veSharded[np], getShardedPair(ve[ib], pm))
			//veSharded[np] = append(veSharded[np], ve[ib])
			ib++
			if ib == lve {
				return
			}
		}
		vNum = ve[ib][0]
		for ib < lve && ve[ib][0] == vNum {
			veSharded[np] = append(veSharded[np], getShardedPair(ve[ib], pm))
			//veSharded[np] = append(veSharded[np], ve[ib])
			ib++
			if ib == lve {
				return
			}
		}
	}
	return
}

type ScalarDissipation struct {
	VtoE              []VertexToElement // Sharded vertex to element map, [2] is [vertID, ElementID_Sharded]
	EtoV              []utils.Matrix    // Sharded Element to Vertex map, Kx3
	Epsilon           []utils.Matrix    // Sharded Np x Kmax, Interpolated from element vertices
	EpsilonScalar     [][]float64       // Sharded scalar value of dissipation, one per element
	DissX, DissY      []utils.Matrix    // Sharded Np x Kmax, Dissipation Flux
	Diss2X, Diss2Y    []utils.Matrix    // Sharded Np x Kmax, Second order X and Y derivative, add together for Divergence
	EpsVertex         []float64         // NVerts x 1, Aggregated (Max) of epsilon surrounding each vertex, Not sharded
	PMap              *PartitionMap     // Partition map for the solution shards in K
	U, UClipped       []utils.Matrix    // Sharded scratch areas for assembly and testing of solution values
	Clipper           utils.Matrix      // Matrix used to clip the topmost mode from the solution polynomial, used in shockfinder
	Element           *DG2D.LagrangeElement2D
	S0, Kappa         float64
	BaryCentricCoords utils.Matrix // A thruple(lam0,lam1,lam2) for interpolation for each interior point, Npx3
}

func (sd *ScalarDissipation) shardEtoV(EtoV utils.Matrix) (ev []utils.Matrix) {
	var (
		pm = sd.PMap
		NP = pm.ParallelDegree
		//KMax, _ = EtoV.Dims()
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

func NewScalarDissipation(kappa float64, dfr *DG2D.DFR2D, pm *PartitionMap) (sd *ScalarDissipation) {
	var (
		el     = dfr.SolutionElement
		NPar   = pm.ParallelDegree
		Np     = el.Np
		order  = float64(el.N)
		NVerts = dfr.VX.Len()
	)
	sd = &ScalarDissipation{
		EpsVertex:     make([]float64, NVerts),
		EpsilonScalar: make([][]float64, NPar),    // Viscosity, constant over the element
		Epsilon:       make([]utils.Matrix, NPar), // Epsilon field, expressed over solution points
		DissX:         make([]utils.Matrix, NPar),
		DissY:         make([]utils.Matrix, NPar),
		Diss2X:        make([]utils.Matrix, NPar),
		Diss2Y:        make([]utils.Matrix, NPar),
		VtoE:          NewVertexToElement(dfr.Tris.EToV).Shard(pm),
		PMap:          pm,
		Element:       el,
		// Sharded working matrices
		U:        make([]utils.Matrix, NPar),
		UClipped: make([]utils.Matrix, NPar),
		S0:       1. / math.Pow(order, 4.),
		Kappa:    0.25,
	}
	sd.EtoV = sd.shardEtoV(dfr.Tris.EToV)
	sd.createInterpolationStencil()
	if kappa != 0. {
		sd.Kappa = kappa
	}
	for np := 0; np < NPar; np++ {
		sd.U[np] = utils.NewMatrix(Np, 1)
		sd.UClipped[np] = utils.NewMatrix(Np, 1)
		Kmax := pm.GetBucketDimension(np)
		sd.Epsilon[np] = utils.NewMatrix(Np, Kmax)
		sd.EpsilonScalar[np] = make([]float64, Kmax)
		sd.DissX[np] = utils.NewMatrix(Np, Kmax)
		sd.DissY[np] = utils.NewMatrix(Np, Kmax)
		sd.Diss2X[np] = utils.NewMatrix(Np, Kmax)
		sd.Diss2Y[np] = utils.NewMatrix(Np, Kmax)
	}
	/*
		The "Clipper" matrix drops the last mode from the polynomial and forms an alternative field of values at the node
		points based on a polynomial with one less term. In other words, if we have a polynomial of degree "p", expressed
		as values at Np node points, multiplying the Node point values vector by Clipper produces an alternative version
		of the node values based on truncating the last polynomial mode.
	*/
	{
		data := make([]float64, Np)
		for i := 0; i < Np; i++ {
			if i != Np-1 {
				data[i] = 1.
			} else {
				data[i] = 0.
			}
		}
		diag := utils.NewDiagMatrix(Np, data)
		sd.Clipper = sd.Element.V.Mul(diag).Mul(sd.Element.Vinv)
	}
	return
}

func (sd *ScalarDissipation) propagateEpsilonMaxToVertices(myThread int, wg *sync.WaitGroup) {
	var (
		VtoE = sd.VtoE[myThread]
		max  = math.Max
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
	wg.Done()
}

type ContinuityLevel uint8

const (
	No ContinuityLevel = iota
	C0
	C1
)

func (sd *ScalarDissipation) AddDissipation(cont ContinuityLevel, myThread int, JdetAll []utils.Matrix, Qall, RHSQall [][4]utils.Matrix) {
	// TODO: Find/Fix memory allocation leak while re-writing this
	var (
		Jdet           = JdetAll[myThread]
		Q              = Qall[myThread]
		DissX, DissY   = sd.DissX[myThread], sd.DissY[myThread]
		Diss2X, Diss2Y = sd.Diss2X[myThread], sd.Diss2Y[myThread]
		RHSQ           = RHSQall[myThread]
		EpsilonScalar  = sd.EpsilonScalar[myThread]
		Epsilon        = sd.Epsilon[myThread]
		Np, KMax       = sd.Element.Np, sd.PMap.GetBucketDimension(myThread)
	)
	for n := 0; n < 4; n++ {
		sd.Element.Dr.Mul(Q[n], DissX) // dQ/dR
		sd.Element.Ds.Mul(Q[n], DissY) // dQ/dS
		switch cont {
		case No:
			for k := 0; k < KMax; k++ {
				for i := 0; i < Np; i++ {
					ind := k + KMax*i
					DissX.DataP[ind] *= EpsilonScalar[k] // Scalar viscosity, constant within each k'th element
					DissY.DataP[ind] *= EpsilonScalar[k]
				}
			}
		case C0:
			// Interpolate epsilon within each element
			sd.linearInterpolateEpsilon(myThread)
			for k := 0; k < KMax; k++ {
				for i := 0; i < Np; i++ {
					ind := k + KMax*i
					DissX.DataP[ind] *= Epsilon.DataP[ind] // C0 viscosity, linear within each k'th element
					DissY.DataP[ind] *= Epsilon.DataP[ind]
				}
			}
		}
		sd.Element.Dr.Mul(DissX, Diss2X)
		sd.Element.Ds.Mul(DissY, Diss2Y)
		Diss2X.Add(Diss2Y) // Compose Divergence in R,S coordinates
		for k := 0; k < KMax; k++ {
			ooJdet2K := 1. / (Jdet.DataP[k] * Jdet.DataP[k]) // Second derivative requires divide by det^2
			for i := 0; i < Np; i++ {
				ind := k + KMax*i
				RHSQ[n].DataP[ind] -= ooJdet2K * Diss2X.DataP[ind]
			}
		}
	}
}

func (sd *ScalarDissipation) linearInterpolateEpsilon(myThread int) {
	var (
		Np, KMax = sd.Element.Np, sd.PMap.GetBucketDimension(myThread)
		Epsilon  = sd.Epsilon[myThread]
		EtoV     = sd.EtoV[myThread]
		R, S     = sd.Element.R, sd.Element.S
	)
	vertLinear := func(r, s float64, f [3]float64) (fi float64) {
		var (
			rLen, sLen     = 2., 2.
			drFrac, dsFrac = (r - (-1.)) / rLen, (s - (-1.)) / sLen
			dfr, dfs       = f[1] - f[0], f[2] - f[0]
		)
		fi = dfr*drFrac + dfs*dsFrac + f[0]
		return
	}
	// Interpolate epsilon within each element
	for k := 0; k < KMax; k++ {
		tri := EtoV.Row(k).DataP
		v := [3]int{int(tri[0]), int(tri[1]), int(tri[2])}
		eps := [3]float64{sd.EpsVertex[v[0]], sd.EpsVertex[v[1]], sd.EpsVertex[v[2]]}
		for i := 0; i < Np; i++ {
			ind := k + KMax*i
			Epsilon.DataP[ind] = vertLinear(R.DataP[i], S.DataP[i], eps)
		}
	}
}

func (sd *ScalarDissipation) baryCentricInterpolateEpsilon(myThread int) {
	var (
		Np, KMax = sd.Element.Np, sd.PMap.GetBucketDimension(myThread)
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
	fld = c.RecombineShardsK(sd.Epsilon)
	return
}

func (sd *ScalarDissipation) GetC0EpsilonPlotField(c *Euler) (fld utils.Matrix) {
	fld = c.RecombineShardsK(sd.Epsilon)
	return
}

func (sd *ScalarDissipation) CalculateElementViscosity(JdetAll []utils.Matrix, Qall [][4]utils.Matrix) {
	var (
		wg = sync.WaitGroup{}
	)
	for np := 0; np < sd.PMap.ParallelDegree; np++ {
		wg.Add(1)
		go func(myThread int) {
			var (
				Rho      = Qall[myThread][0]
				Jdet     = JdetAll[myThread]
				Eps      = sd.EpsilonScalar[myThread]
				Kmax     = sd.PMap.GetBucketDimension(myThread)
				U        = sd.U[myThread]
				UClipped = sd.UClipped[myThread]
			)
			for k := 0; k < Kmax; k++ {
				var (
					eps0        = math.Sqrt(2.*Jdet.DataP[k]) / float64(sd.Element.N)
					Se          = math.Log10(sd.moment(k, Kmax, U, UClipped, Rho))
					left, right = sd.S0 - sd.Kappa, sd.S0 + sd.Kappa
					oo2kappa    = 0.5 / sd.Kappa
				)
				switch {
				case Se < left:
					Eps[k] = 0.
				case Se >= left && Se < right:
					Eps[k] = 0.5 * eps0 * (1. + math.Sin(math.Pi*oo2kappa*(Se-sd.S0)))
				case Se >= right:
					Eps[k] = eps0
				}
			}
			wg.Done()
		}(np)
	}
	wg.Wait()
	for np := 0; np < sd.PMap.ParallelDegree; np++ {
		wg.Add(1)
		go sd.propagateEpsilonMaxToVertices(np, &wg)
	}
	wg.Wait()
}

func (sd *ScalarDissipation) moment(k, Kmax int, U, UClipped, Rho utils.Matrix) (m float64) {
	var (
		Np            = sd.Element.Np
		UD, UClippedD = U.DataP, UClipped.DataP
	)
	for i := 0; i < Np; i++ {
		ind := k + i*Kmax
		U.DataP[i] = Rho.DataP[ind]
	}
	/*
		Evaluate the L2 moment of (q - qalt) over the element, where qalt is the truncated version of q
		Here we don't bother using quadrature, we do a simple sum
	*/
	UClipped = sd.Clipper.Mul(U, UClipped)
	for i := 0; i < Np; i++ {
		t1 := UD[i] - UClippedD[i]
		m += t1 * t1 / (UD[i] * UD[i])
	}
	return
}

func (sd *ScalarDissipation) createInterpolationStencil() {
	var (
		Np    = sd.Element.Np
		R, S  = sd.Element.R, sd.Element.S
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
			//sd.BaryCentricCoords.Set(i, ii, LAM.DataP[ii])
		}
	}
}

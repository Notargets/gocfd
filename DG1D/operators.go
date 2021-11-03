package DG1D

import (
	"math"

	"github.com/notargets/gocfd/utils"
)

func (el Elements1D) SlopeLimitN(U utils.Matrix, M float64) (ULim utils.Matrix) {
	var (
		Uh    = el.Vinv.Mul(U)
		eps0  = 1.0e-8
		nr, _ = U.Dims()
	)
	Uh.SetRange(1, -1, 0, -1, 0)
	Uh = el.V.Mul(Uh)
	vk := Uh.Row(0)

	// End values of each element
	ue1 := U.Row(0)
	ue2 := U.Row(nr - 1)

	// Cell averages
	vkm1 := vk.Subset(0, 0).Concat(vk.Subset(0, -2))
	vkp1 := vk.Subset(1, -1).Concat(vk.Subset(-1, -1))

	// Apply reconstruction to find elements in need of limiting
	vm1 := vk.Copy().Subtract(vkm1)
	vp1 := vkp1.Copy().Subtract(vk)
	var ve1, ve2 utils.Vector
	if M == 0 {
		ve1 = vk.Copy().Subtract(Minmod(vk.Copy().Subtract(ue1), vm1, vp1))
		ve2 = vk.Copy().Add(Minmod(ue2.Copy().Subtract(vk), vm1, vp1))
	} else {
		h := el.X.Row(0).AtVec(1) - el.X.Row(0).AtVec(0)
		ve1 = vk.Copy().Subtract(MinmodB(M, h, vk.Copy().Subtract(ue1), vm1, vp1))
		ve2 = vk.Copy().Add(MinmodB(M, h, ue2.Copy().Subtract(vk), vm1, vp1))
	}

	ids := ve1.Subtract(ue1).FindOr(utils.Greater, eps0, true, ve2.Subtract(ue2))

	ULim = U.Copy()
	if ids.Len() != 0 {
		//fmt.Printf("ids = %v\n", ids.ToIndex())
		idsI := ids.ToIndex()
		// We need to limit the elements in the index
		// Create a piecewise linear solution for limiting the elements in the index
		uhl := el.Vinv.Mul(U.SliceCols(idsI))
		uhl.SetRange(2, -1, 0, -1, 0) // Set all polynomial coefficients higher than linear to 0
		ul := el.V.Mul(uhl)
		// Apply slope limiter to specified elements
		ULim.AssignColumns(idsI, el.SlopeLimitLin(ul, el.X.SliceCols(idsI), vkm1.SubsetIndex(idsI), vk.SubsetIndex(idsI), vkp1.SubsetIndex(idsI)))
	}
	return
}

func (el Elements1D) SlopeLimitLin(ul, xl utils.Matrix, vm1, v0, vp1 utils.Vector) (ULim utils.Matrix) {
	var (
		Np       = el.Np
		ones     = utils.NewVectorConstant(Np, 1)
		h        = xl.Row(Np - 1).Subtract(xl.Row(0))
		x0       = ones.Outer(xl.Row(0).Add(h.Copy().Scale(0.5)))
		hNScaled = ones.Outer(h).POW(-1).Scale(2)
		ux       = hNScaled.ElMul(el.Dr.Copy().Mul(ul))
	)
	ULim = ones.Outer(v0).Add(xl.Subtract(x0).ElMul(ones.Outer(Minmod(ux.Row(0), vp1.Subtract(v0).ElDiv(h), v0.Subtract(vm1).ElDiv(h)))))
	return
}

func Minmod(vecs ...utils.Vector) (R utils.Vector) {
	/*
		Computes minmod across a group of vectors
		    Input: Ainv, B, C, length N
				For each element in Ainv, B, C, compose a vector like {a1, b1, c1} and set r1 = minmod(a1,b1,c1)
			Output: R, length N
	*/
	var (
		W     = len(vecs)
		dataV = make([]float64, W)
		N     = vecs[0].Len()
		dataR = make([]float64, N)
	)
	for i := 0; i < N; i++ {
		for j := 0; j < W; j++ {
			dataV[j] = vecs[j].AtVec(i)
		}
		dataR[i] = minmod(dataV)
	}
	R = utils.NewVector(N, dataR)
	return
}

func minmod(a []float64) (r float64) {
	var (
		rMin = math.Abs(a[0])
	)
	var signSum int
	for _, val := range a {
		if math.Signbit(val) {
			signSum -= 1
		} else {
			signSum += 1
		}
		rMin = math.Min(rMin, math.Abs(val))
	}
	ss := int(math.Abs(float64(signSum)))
	sign := signSum / ss
	switch ss {
	case len(a):
		return float64(sign) * rMin
	default:
		return 0
	}
}

func MinmodB(M, h float64, vecs ...utils.Vector) (R utils.Vector) {
	/*
		Computes minmodB across a group of vectors
		    Input: Ainv, B, C, length N
				For each element in Ainv, B, C, compose a vector like {a1, b1, c1} and set r1 = minmod(a1,b1,c1)
			Output: R, length N
	*/
	var (
		W     = len(vecs)
		dataV = make([]float64, W)
		N     = vecs[0].Len()
		dataR = make([]float64, N)
	)
	for i := 0; i < N; i++ {
		dataR[i] = vecs[0].RawVector().Data[i]
	}
	// Check for values higher than our limit in the first vector
	factor := M * utils.POW(h, 2)
	idsV := vecs[0].Find(utils.Greater, factor, true)
	if idsV.Len() != 0 {
		ids := idsV.ToIndex()
		for _, i := range ids {
			for j := 0; j < W; j++ {
				dataV[j] = vecs[j].AtVec(i)
			}
			dataR[i] = minmod(dataV)
		}
	}
	R = utils.NewVector(N, dataR)
	return
}

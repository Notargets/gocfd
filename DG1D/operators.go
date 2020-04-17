package DG1D

import (
	"math"

	"github.com/notargets/gocfd/utils"
)

func (el Elements1D) SlopeLimitN(U utils.Matrix) (ULim utils.Matrix) {
	var (
		Uh   = el.Vinv.Mul(U)
		eps0 = 1.0e-8
	)
	Uh.SetRange(1, -1, 0, -1, 0)
	Uh = el.V.Mul(Uh)
	vk := Uh.Row(0)

	// End values of each element
	ue1 := U.Row(0)
	ue2 := U.Row(-1)

	// Cell averages
	vkm1 := vk.Subset(0, 1).Concat(vk.Subset(1, -2))
	vkp1 := vk.Subset(1, -1).Concat(vk.Subset(-1, -1))

	// Apply reconstruction to find elements in need of limiting
	vm1 := vk.Copy().Subtract(vkm1)
	vp1 := vkp1.Copy().Subtract(vk)
	ve1 := vk.Copy().Subtract(Minmod(vk.Copy().Subtract(ue1), vm1, vp1))
	ve2 := vk.Copy().Add(Minmod(ue2.Copy().Subtract(vk), vm1, vp1))
	ids := ve1.Subtract(ue1).Find(utils.Greater, eps0, true)
	ids = ids.Concat(ve2.Subtract(ue2).Find(utils.Greater, eps0, true))
	return
}

func Minmod(vecs ...utils.Vector) (R utils.Vector) {
	/*
		Computes minmod across a group of vectors
		    Input: A, B, C, length N
				For each element in A, B, C, compose a vector like {a1, b1, c1} and set r1 = minmod(a1,b1,c1)
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
		rMin = a[0]
	)
	for _, val := range a {
		rMin = math.Min(rMin, val)
	}
	return math.Max(0, rMin)
}

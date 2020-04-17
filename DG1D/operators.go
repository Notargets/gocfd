package DG1D

import "github.com/notargets/gocfd/utils"

func (el Elements1D) SlopeLimitN(U utils.Matrix) (ULim utils.Matrix) {
	var (
		Uh = el.Vinv.Mul(U)
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
	_, _, _, _ = ue1, ue2, vkm1, vkp1
	return
}

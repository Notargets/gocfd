package Euler2D

import "github.com/notargets/gocfd/utils"

func Get4DP(Q [4]utils.Matrix) (qD [4][]float64) {
	qD = [4][]float64{Q[0].Data(), Q[1].Data(), Q[2].Data(), Q[3].Data()}
	return
}

func (c *Euler) split1D(iMax, threadNum int) (istart, iend int) {
	var (
		Npart = iMax / c.ParallelDegree
	)
	istart = threadNum * Npart
	iend = istart + Npart
	if threadNum == c.ParallelDegree-1 {
		iend = iMax
	}
	return
}

func (c *Euler) GetQQ(k, i int, Q [4][]float64) (qq [4]float64) {
	var (
		Kmax = c.dfr.K
		ind  = k + Kmax*i
	)
	qq = [4]float64{Q[0][ind], Q[1][ind], Q[2][ind], Q[3][ind]}
	return
}

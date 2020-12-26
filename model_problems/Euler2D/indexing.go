package Euler2D

import (
	"github.com/notargets/gocfd/utils"
)

func Get4DP(Q [4]utils.Matrix) (qD [4][]float64) {
	qD = [4][]float64{Q[0].Data(), Q[1].Data(), Q[2].Data(), Q[3].Data()}
	return
}

func (c *Euler) split1D(iMax, threadNum int) (istart, iend int) {
	// This routine splits one dimension into c.ParallelDegree pieces, with a maximum imbalance of one item
	var (
		Npart            = iMax / (c.ParallelDegree)
		startAdd, endAdd int
		remainder        int
	)
	remainder = iMax % c.ParallelDegree
	if remainder != 0 { // spread the remainder over the first chunks evenly
		if threadNum+1 > remainder {
			startAdd = remainder
			endAdd = 0
		} else {
			startAdd = threadNum
			endAdd = 1
		}
	}
	istart = threadNum*Npart + startAdd
	iend = istart + Npart + endAdd
	return
}

func (c *Euler) GetKSplitRange(threadnum int) (kStart, kEnd int) {
	// Returns the global k range for thread number "threadnum"
	kStart, kEnd = c.split1D(c.dfr.K, threadnum)
	return
}
func (c *Euler) GetKSplitOffset(threadnum int) (kOffset int) {
	kOffset, _ = c.GetKSplitRange(threadnum)
	return
}
func (c *Euler) GetKSplitMaxK(threadnum int) (kMax int) {
	var (
		k1, k2 = c.GetKSplitRange(threadnum)
	)
	kMax = k2 - k1
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

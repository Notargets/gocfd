package Euler2D

import (
	"runtime"

	"github.com/notargets/gocfd/utils"
)

func (c *Euler) ShardByK(A utils.Matrix) (pA []utils.Matrix) {
	var (
		NP      = c.Partitions.ParallelDegree
		Imax, _ = A.Dims()
		aD      = A.Data()
	)
	pA = make([]utils.Matrix, NP)
	for np := 0; np < NP; np++ {
		kMin, kMax := c.Partitions.GetBucketRange(np)
		Kmax := c.Partitions.GetBucketDimension(np)
		pA[np] = utils.NewMatrix(Imax, Kmax)
		pAD := pA[np].Data()
		for k := kMin; k < kMax; k++ {
			pk := k - kMin
			for i := 0; i < Imax; i++ {
				ind := k + i*c.dfr.K
				pind := pk + Kmax*i
				pAD[pind] = aD[ind]
			}
		}
	}
	return
}

func (c *Euler) ShardByKTranspose(A utils.Matrix) (pA []utils.Matrix) {
	var (
		NP      = c.Partitions.ParallelDegree
		_, Imax = A.Dims()
		aD      = A.Data()
	)
	pA = make([]utils.Matrix, NP)
	for np := 0; np < NP; np++ {
		kMin, kMax := c.Partitions.GetBucketRange(np)
		Kmax := c.Partitions.GetBucketDimension(np)
		pA[np] = utils.NewMatrix(Kmax, Imax)
		pAD := pA[np].Data()
		for k := kMin; k < kMax; k++ {
			pk := k - kMin
			for i := 0; i < Imax; i++ {
				//ind := k + i*c.dfr.K
				//pind := pk + Kmax*i
				ind := i + k*Imax
				pind := i + pk*Imax
				pAD[pind] = aD[ind]
			}
		}
	}
	return
}

func (c *Euler) RecombineShardsK(pA []utils.Matrix) (A utils.Matrix) {
	var (
		NP      = c.Partitions.ParallelDegree
		_, Imax = pA[0].Dims()
	)
	A = utils.NewMatrix(Imax, c.dfr.K)
	aD := A.Data()
	for np := 0; np < NP; np++ {
		kMin, kMax := c.Partitions.GetBucketRange(np)
		Kmax := c.Partitions.GetBucketDimension(np)
		pAD := pA[np].Data()
		for k := kMin; k < kMax; k++ {
			pk := k - kMin
			for i := 0; i < Imax; i++ {
				ind := k + i*c.dfr.K
				pind := pk + Kmax*i
				aD[ind] = pAD[pind]
			}
		}
	}
	return
}

func (c *Euler) RecombineShardsKBy4(pA [][4]utils.Matrix) (A [4]utils.Matrix) {
	var (
		NP      = c.Partitions.ParallelDegree
		Imax, _ = pA[0][0].Dims()
	)
	for n := 0; n < 4; n++ {
		A[n] = utils.NewMatrix(Imax, c.dfr.K)
		aD := A[n].Data()
		for np := 0; np < NP; np++ {
			kMin, kMax := c.Partitions.GetBucketRange(np)
			Kmax := c.Partitions.GetBucketDimension(np)
			pAD := pA[np][n].Data()
			for k := kMin; k < kMax; k++ {
				pk := k - kMin
				for i := 0; i < Imax; i++ {
					ind := k + i*c.dfr.K
					pind := pk + Kmax*i
					aD[ind] = pAD[pind]
				}
			}
		}
	}
	return
}

func (c *Euler) SetParallelDegree(ProcLimit, Kmax int) {
	if ProcLimit != 0 {
		c.ParallelDegree = ProcLimit
	} else {
		c.ParallelDegree = runtime.NumCPU()
	}
	runtime.GOMAXPROCS(runtime.NumCPU())
	if c.ParallelDegree > Kmax {
		c.ParallelDegree = 1
	}
	c.Partitions = NewPartitionMap(c.ParallelDegree, Kmax)
}

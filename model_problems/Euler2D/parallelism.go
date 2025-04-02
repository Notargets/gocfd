package Euler2D

import (
	"fmt"
	"math"
	"runtime"

	"github.com/notargets/gocfd/utils"
)

func (c *Euler) ShardByK(A utils.Matrix) (pA []utils.Matrix) {
	var (
		NP      = c.Partitions.ParallelDegree
		Imax, _ = A.Dims()
		aD      = A.DataP
	)
	pA = make([]utils.Matrix, NP)
	for np := 0; np < NP; np++ {
		kMin, kMax := c.Partitions.GetBucketRange(np)
		Kmax := c.Partitions.GetBucketDimension(np)
		pA[np] = utils.NewMatrix(Imax, Kmax)
		pAD := pA[np].DataP
		for k := kMin; k < kMax; k++ {
			pk := k - kMin
			for i := 0; i < Imax; i++ {
				ind := k + i*c.DFR.K
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
		aD      = A.DataP
	)
	pA = make([]utils.Matrix, NP)
	for np := 0; np < NP; np++ {
		kMin, kMax := c.Partitions.GetBucketRange(np)
		Kmax := c.Partitions.GetBucketDimension(np)
		pA[np] = utils.NewMatrix(Kmax, Imax)
		pAD := pA[np].DataP
		for k := kMin; k < kMax; k++ {
			pk := k - kMin
			for i := 0; i < Imax; i++ {
				// ind := k + i*c.DFR.K
				// pind := pk + Kmax*i
				ind := i + k*Imax
				pind := i + pk*Imax
				pAD[pind] = aD[ind]
			}
		}
	}
	return
}

func (c *Euler) RecombineShardsK(pA []utils.Matrix, A *utils.Matrix) {
	var (
		NP      = c.Partitions.ParallelDegree
		Imax, _ = pA[0].Dims()
	)
	if A.IsEmpty() {
		*A = utils.NewMatrix(Imax, c.DFR.K)
	}
	aD := A.DataP
	for np := 0; np < NP; np++ {
		kMin, kMax := c.Partitions.GetBucketRange(np)
		Kmax := c.Partitions.GetBucketDimension(np)
		pAD := pA[np].DataP
		for k := kMin; k < kMax; k++ {
			pk := k - kMin
			for i := 0; i < Imax; i++ {
				ind := k + i*c.DFR.K
				pind := pk + Kmax*i
				aD[ind] = pAD[pind]
			}
		}
	}
	return
}

func (c *Euler) RecombineShardsKBy4(pA [][4]utils.Matrix, A *[4]utils.Matrix) {
	var (
		NP = c.Partitions.ParallelDegree
	)
	ppA := make([]utils.Matrix, NP)
	for n := 0; n < 4; n++ {
		for np := 0; np < NP; np++ {
			ppA[np] = pA[np][n]
		}
		c.RecombineShardsK(ppA, &A[n])
	}
	return
}

func (c *Euler) SetParallelDegree(ProcLimit, Kmax int) {
	/*
		Here we set the number of partitions we'll use for parallelism.
		Numerical experiments suggest that there should be a maximum of 185 elements per partition,
		so we'll use that as our baseline
	*/
	var (
		ParallelDegree int
	)
	if ProcLimit != 0 {
		ParallelDegree = ProcLimit
	} else {
		ParallelDegree = Kmax / 185
		ParallelDegree = int(math.Max(float64(runtime.NumCPU()), float64(ParallelDegree)))
	}
	runtime.GOMAXPROCS(runtime.NumCPU())
	if ParallelDegree > Kmax {
		ParallelDegree = 1
	}
	c.Partitions = NewPartitionMap(ParallelDegree, Kmax)
}

func (c *Euler) PartitionEdges() {
	var (
		NPar                               = c.Partitions.ParallelDegree
		SharedEdges, BCEdges, PhantomEdges EdgeKeySlice
	)
	// First, separate edges into 3 groups
	for en, e := range c.DFR.Tris.Edges {
		switch e.NumConnectedTris {
		case 0:
			PhantomEdges = append(PhantomEdges, en)
		case 1:
			BCEdges = append(BCEdges, en)
		case 2:
			SharedEdges = append(SharedEdges, en)
		}
	}
	if len(SharedEdges) == 0 && len(BCEdges) == 0 {
		err := fmt.Errorf("Number of edges should be > 0, have Shared[%d], BC[%d], Phantom[%d]\n",
			len(SharedEdges), len(BCEdges), len(PhantomEdges))
		panic(err)
	}
	c.SortedEdgeKeys = make([]EdgeKeySlice, NPar)
	pmS := NewPartitionMap(NPar, len(SharedEdges))
	pmB := NewPartitionMap(NPar, len(BCEdges))
	pmP := NewPartitionMap(NPar, len(PhantomEdges))
	for np := 0; np < NPar; np++ {
		SSize := pmS.GetBucketDimension(np)
		BSize := pmB.GetBucketDimension(np)
		PSize := pmP.GetBucketDimension(np)
		c.SortedEdgeKeys[np] = make(EdgeKeySlice, SSize+BSize+PSize)
	}
	for np := 0; np < NPar; np++ {
		var ii int
		if len(PhantomEdges) != 0 {
			PMin, PMax := pmP.GetBucketRange(np)
			for i := PMin; i < PMax; i++ {
				c.SortedEdgeKeys[np][ii] = PhantomEdges[i]
				ii++
			}
		}
		if len(BCEdges) != 0 {
			BMin, BMax := pmB.GetBucketRange(np)
			for i := BMin; i < BMax; i++ {
				c.SortedEdgeKeys[np][ii] = BCEdges[i]
				ii++
			}
		}
		if len(SharedEdges) != 0 {
			SMin, SMax := pmS.GetBucketRange(np)
			for i := SMin; i < SMax; i++ {
				c.SortedEdgeKeys[np][ii] = SharedEdges[i]
				ii++
			}
		}
	}
}

func (c *Euler) PartitionEdgesByK() {
	var (
		pm   = c.Partitions
		NPar = pm.ParallelDegree
	)
	c.SortedEdgeKeys = make([]EdgeKeySlice, NPar)
	// First, separate edges into 3 groups
	for en, e := range c.DFR.Tris.Edges {
		bn, _, _ := pm.GetBucket(int(e.ConnectedTris[0]))
		c.SortedEdgeKeys[bn] = append(c.SortedEdgeKeys[bn], en)
	}
}

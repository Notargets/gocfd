package Euler2D

type PartitionMap struct {
	MaxIndex       int // MaxIndex is partitioned into ParallelDegree partitions
	ParallelDegree int
	Partitions     [][2]int // Beginning and end index of partitions
}

func NewPartitionMap(ParallelDegree, maxIndex int) (pm *PartitionMap) {
	pm = &PartitionMap{
		MaxIndex:       maxIndex,
		ParallelDegree: ParallelDegree,
		Partitions:     make([][2]int, ParallelDegree),
	}
	for n := 0; n < ParallelDegree; n++ {
		pm.Partitions[n] = pm.Split1D(n)
	}
	return
}

func (pm *PartitionMap) GetBucket(kDim int) (bucketNum, min, max int) {
	_, bucketNum, min, max = pm.getBucketWithTryCount(kDim)
	return
}

func (pm *PartitionMap) getBucketWithTryCount(kDim int) (tryCount, bucketNum, min, max int) {
	// Initial guess
	bucketNum = int(float64(pm.ParallelDegree*kDim) / float64(pm.MaxIndex))
	for !(pm.Partitions[bucketNum][0] <= kDim && pm.Partitions[bucketNum][1] > kDim) {
		if pm.Partitions[bucketNum][0] > kDim {
			bucketNum--
		} else {
			bucketNum++
		}
		if bucketNum == -1 || bucketNum == pm.ParallelDegree {
			return 0, -1, 0, 0
		}
		tryCount++
	}
	min, max = pm.Partitions[bucketNum][0], pm.Partitions[bucketNum][1]
	/*
		if tryCount != 0 {
			fmt.Printf("bn, kDim, maxIndex, ParallelDegree, tryCount = %d, %d, %d, %d, %d\n",
				bucketNum, kDim, pm.MaxIndex, pm.ParallelDegree, tryCount)
		}
	*/
	return
}

func (pm *PartitionMap) GetBucketRange(bucketNum int) (kMin, kMax int) {
	kMin, kMax = pm.Partitions[bucketNum][0], pm.Partitions[bucketNum][1]
	return
}

func (pm *PartitionMap) GetLocalK(baseK int) (k, Kmax, bn int) {
	var (
		kmin, kmax int
	)
	bn, kmin, kmax = pm.GetBucket(baseK)
	Kmax = kmax - kmin
	k = baseK - kmin
	return
}

func (pm *PartitionMap) GetBucketDimension(bucketNum int) (kMax int) {
	var (
		k1, k2 = pm.GetBucketRange(bucketNum)
	)
	kMax = k2 - k1
	return
}

func (pm *PartitionMap) Split1D(threadNum int) (bucket [2]int) {
	// This routine splits one dimension into c.ParallelDegree pieces, with a maximum imbalance of one item
	var (
		Npart            = pm.MaxIndex / (pm.ParallelDegree)
		startAdd, endAdd int
		remainder        int
	)
	remainder = pm.MaxIndex % pm.ParallelDegree
	if remainder != 0 { // spread the remainder over the first chunks evenly
		if threadNum+1 > remainder {
			startAdd = remainder
			endAdd = 0
		} else {
			startAdd = threadNum
			endAdd = 1
		}
	}
	bucket[0] = threadNum*Npart + startAdd
	bucket[1] = bucket[0] + Npart + endAdd
	return
}

func (c *Euler) GetQQ(k, Kmax, i int, Q [4][]float64) (qq [4]float64) {
	var (
		ind = k + Kmax*i
	)
	qq = [4]float64{Q[0][ind], Q[1][ind], Q[2][ind], Q[3][ind]}
	return
}

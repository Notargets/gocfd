package utils

import (
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestEuler_Indexing(t *testing.T) {
	{ // Test PartitionMap
		getHisto := func(K, Np int) (histo map[int]int) {
			pm := NewPartitionMap(Np, K)
			histo = make(map[int]int)
			for np := 0; np < pm.ParallelDegree; np++ {
				maxK := pm.GetBucketDimension(np)
				histo[maxK]++
			}
			return
		}
		getTotal := func(histo map[int]int) (total int) {
			for key, count := range histo {
				total += key * count
			}
			return
		}
		assert.Equal(t, map[int]int{0: 30, 1: 2}, getHisto(2, 32))
		assert.Equal(t, map[int]int{1: 32}, getHisto(32, 32))
		assert.Equal(t, map[int]int{8: 32}, getHisto(256, 32))
		assert.Equal(t, map[int]int{8: 1, 9: 31}, getHisto(287, 32))
		assert.Equal(t, 287, getTotal(getHisto(287, 32)))
		for n := 64; n < 10000; n++ {
			// for n := 64; n < 10000; n++ {
			// n := 64
			// {
			var (
				keys   [2]float64
				keyNum int
			)
			histo := getHisto(n, 32)
			for key := range histo {
				keys[keyNum] = float64(key)
				keyNum++
			}
			if keyNum == 2 {
				assert.Equal(t, 1., math.Abs(keys[0]-keys[1])) // Maximum imbalance of 1
			}
			// fmt.Printf("keys = %v, histo[%d] = %v\n", keys, n, histo)
			assert.Equal(t, n, getTotal(histo))
		}
	}
	{ // Test inverted bucket probe - find bucket that contains index (efficiently)
		for maxIndex := 10; maxIndex < 1000; maxIndex++ {
			pm := NewPartitionMap(5, maxIndex)
			for k := 0; k < maxIndex; k++ {
				tryCount, bn, min, max := pm.getBucketWithTryCount(k)
				mmin, mmax := pm.GetBucketRange(bn)
				assert.True(t, k >= min && k < max && min == mmin && max == mmax && tryCount <= 1)
			}
		}
	}
}

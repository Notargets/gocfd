package Euler2D

import (
	"math"
	"testing"

	"github.com/stretchr/testify/assert"

	"github.com/notargets/gocfd/DG2D"
)

func TestEuler_Indexing(t *testing.T) {
	{
		getHisto := func(K, Np int) (histo map[int]int) {
			histo = make(map[int]int)
			d := &DG2D.DFR2D{K: K}
			c := &Euler{dfr: d, ParallelDegree: Np}
			for np := 0; np < c.ParallelDegree; np++ {
				maxK := c.GetKSplitMaxK(np)
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
			//for n := 64; n < 10000; n++ {
			//n := 64
			//{
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
			//fmt.Printf("keys = %v, histo[%d] = %v\n", keys, n, histo)
			assert.Equal(t, n, getTotal(histo))
		}
		//fmt.Printf("histogram = %v\n", getHisto(287, 32))
	}
}

package types

import (
	"fmt"
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestTypes(t *testing.T) {
	{ // Test packed int for edge labeling
		en := NewEdgeKey([2]int{1, 0})
		assert.Equal(t, EdgeKey(1<<32), en)
		assert.Equal(t, [2]int{0, 1}, en.GetVertices(false))

		en = NewEdgeKey([2]int{0, 1})
		assert.Equal(t, EdgeKey(1<<32), en)
		assert.Equal(t, [2]int{0, 1}, en.GetVertices(false))

		en = NewEdgeKey([2]int{0, 10})
		assert.Equal(t, EdgeKey(10*(1<<32)), en)
		assert.Equal(t, [2]int{0, 10}, en.GetVertices(false))

		en = NewEdgeKey([2]int{100, 0})
		assert.Equal(t, EdgeKey(100*(1<<32)), en)
		assert.Equal(t, [2]int{0, 100}, en.GetVertices(false))

		en = NewEdgeKey([2]int{100, 1})
		assert.Equal(t, EdgeKey(100*(1<<32)+1), en)
		assert.Equal(t, [2]int{1, 100}, en.GetVertices(false))

		en = NewEdgeKey([2]int{100, 100001})
		assert.Equal(t, EdgeKey(100001*(1<<32)+100), en)
		assert.Equal(t, [2]int{100, 100001}, en.GetVertices(false))

		// Test maximum/minimum indices
		en = NewEdgeKey([2]int{1, 1<<32 - 1})
		assert.Equal(t, EdgeKey((1<<32-1)<<32+1), en)
		assert.Equal(t, [2]int{1, 1<<32 - 1}, en.GetVertices(false))

		en = NewEdgeKey([2]int{1<<32 - 1, 1<<32 - 1})
		assert.Equal(t, EdgeKey(1<<64-1), en)
		assert.Equal(t, [2]int{1<<32 - 1, 1<<32 - 1}, en.GetVertices(false))

		en = NewEdgeKey([2]int{1<<32 - 1, 1})
		assert.Equal(t, EdgeKey((1<<32-1)<<32+1), en)
		assert.Equal(t, [2]int{1, 1<<32 - 1}, en.GetVertices(false))
	}
	{
		tokens := []string{"WALL", "Periodic-1", "Periodic-2", "Wall-22", "Wall-top", "Neuman-10"}
		flags := []BCFLAG{BC_Wall, BC_Periodic, BC_Periodic, BC_Wall, BC_Wall, BC_Neuman}
		labels := []string{"", "1", "2", "22", "top", "10"}
		for i, token := range tokens {
			bt := NewBCTAG(token)
			fmt.Printf("bt = %s, bcflag = %v\n", bt, bt.GetFLAG().String())
			assert.Equal(t, flags[i], bt.GetFLAG())
			assert.Equal(t, labels[i], bt.GetLabel())
		}
	}
	{ // Test curve ordering/reordering
		// combo of int1, int2 equals edges
		//ind1 := []int{40, 39, 38, 37}
		//ind2 := []int{0, 40, 39, 38}
		ind1 := []int{38, 40, 39, 37}
		ind2 := []int{39, 0, 40, 38}
		edges := make(Curve, len(ind1))
		for i := range ind1 {
			edges[i] = NewEdgeInt([2]int{ind1[i], ind2[i]})
		}
		//fmt.Printf("before reordering = \n")
		//edges.Print()
		var unordered bool
		edges, unordered = edges.ReOrder(false)
		assert.Equal(t, true, unordered)
		assert.Equal(t, [2]int{40, 0}, edges[0].GetVertices())
		assert.Equal(t, [2]int{37, 38}, edges[3].GetVertices())
		//fmt.Printf("reordered edges = \n")
		//edges.Print()
		edges, unordered = edges.ReOrder(true)
		assert.Equal(t, false, unordered)
		//fmt.Printf("reverse reordered edges = \n")
		//edges.Print()
		assert.Equal(t, [2]int{37, 38}, edges[0].GetVertices())
		assert.Equal(t, [2]int{40, 0}, edges[3].GetVertices())
	}
}

package types

import (
	"fmt"
	"math"
)

type EdgeNumber uint64

func NewEdgeNumber(verts [2]int) (packed EdgeNumber) {
	// This packs two index coordinates into two 32 bit unsigned integers to act as a hash and an indirect access method
	var (
		limit = math.MaxUint32
	)
	for _, vert := range verts {
		if vert < 0 || vert > limit {
			panic(fmt.Errorf("unable to pack two ints into a uint64, have %d and %d as inputs",
				verts[0], verts[1]))
		}
	}
	var i1, i2 int
	if verts[0] < verts[1] {
		i1, i2 = verts[0], verts[1]
	} else {
		i1, i2 = verts[1], verts[0]
	}
	packed = EdgeNumber(i1 + i2<<32)
	return
}

func (en EdgeNumber) GetVertices(rev bool) (verts [2]int) {
	var (
		enTmp EdgeNumber
	)
	enTmp = en >> 32
	verts[1] = int(enTmp)
	verts[0] = int(en - enTmp*(1<<32))
	if rev {
		verts[0], verts[1] = verts[1], verts[0]
	}
	return
}

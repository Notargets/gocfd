package types

import (
	"fmt"
	"math"
)

/*
EdgeNumber is an always positive number that stores an edge's vertices as indices in a way that can be compared
An edge between vertices [4] and [0] will always be stored as [0,4], in the ascending order of the index values
*/
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
	if verts[0] <= verts[1] {
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

/*
An Edge stores the edge vertices in the original order of the vertices, so that it can be recovered with it's direction
*/
type EdgeInt int64

func NewEdgeInt(verts [2]int) (packed EdgeInt) {
	// This packs two index coordinates into two 31 bit unsigned integers to act as a hash and an indirect access method
	var (
		limit = math.MaxUint32 >> 1 // leaves room for the sign bit of an int64
		sign  bool
	)
	for _, vert := range verts {
		if vert < 0 || vert > limit {
			panic(fmt.Errorf("unable to pack two ints into a uint64, have %d and %d as inputs",
				verts[0], verts[1]))
		}
	}
	var i1, i2 int
	if verts[0] <= verts[1] {
		i1, i2 = verts[0], verts[1]
	} else {
		sign = true
		i1, i2 = verts[1], verts[0]
	}
	packed = EdgeInt(i1 + i2<<32)
	if sign {
		packed = -packed
	}
	return
}

func (e EdgeInt) GetVertices() (verts [2]int) {
	var (
		eTmp EdgeInt
		sign bool
	)
	if e < 0 {
		sign = true
		e = -e
	}
	eTmp = e >> 32
	verts[1] = int(eTmp)
	verts[0] = int(e - eTmp*(1<<32))
	if sign {
		verts[0], verts[1] = verts[1], verts[0]
	}
	return
}

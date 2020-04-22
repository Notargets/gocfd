package utils

import (
	"strconv"
	"strings"
)

type R1 struct {
	Max int
}

func NewR1(imax int) R1 {
	return R1{imax}
}

func (r R1) Range(dimI interface{}) Index {
	var (
		i1, i2 = ParseDim(dimI, r.Max)
	)
	return NewRange(i1, i2)
}

type R2 struct {
	Ir, Jr R1
}

func NewR2(imax, jmax int) R2 {
	return R2{
		NewR1(imax),
		NewR1(jmax),
	}
}

func (r R2) Range(dimI, dimJ interface{}) (I Index) {
	var (
		i1, i2 = ParseDim(dimI, r.Ir.Max)
		j1, j2 = ParseDim(dimJ, r.Jr.Max)
		_, nj  = r.Ir.Max, r.Jr.Max
	)
	size := (i2 - i1) * (j2 - j1)
	I = NewIndex(size)
	var ind int
	for j := j1; j < j2; j++ {
		for i := i1; i < i2; i++ {
			I[ind] = j + nj*i // Column Major
			ind++
		}
	}
	return
}

type R3 struct {
	Ir, Jr, Kr R1
}

func NewR3(imax, jmax, kmax int) R3 {
	return R3{
		NewR1(imax),
		NewR1(jmax),
		NewR1(kmax),
	}
}

func (r R3) Range(dimI, dimJ, dimK interface{}) (I Index) {
	var (
		i1, i2    = ParseDim(dimI, r.Ir.Max)
		j1, j2    = ParseDim(dimJ, r.Jr.Max)
		k1, k2    = ParseDim(dimK, r.Kr.Max)
		ni, nj, _ = r.Ir.Max, r.Jr.Max, r.Kr.Max
	)
	size := (i2 - i1) * (j2 - j1) * (k2 - k1)
	I = NewIndex(size)
	var ind int
	for k := k1; k < k2; k++ {
		for j := j1; j < j2; j++ {
			for i := i1; i < i2; i++ {
				I[ind] = j + nj*(i+ni*k) // Column Major
				ind++
			}
		}
	}
	return
}

func ParseDim(dimI interface{}, max int) (i1, i2 int) {
	/*
		Converts phrases including:
			":"   = full range, from 0 to max (loop indexing)
			"end" = last index, from max-1, max
			"N"   = middle index, from N-1, N
		   	N     = middle index, from N-1, N
		    "2:N" = range, from 2 to N (loop indexing)
		   	":N"  = range, from 0 to N (loop indexing)
		   	"N:"  = range, from N to max-1 (loop indexing)
	*/
	switch dim := dimI.(type) {
	case string:
		switch dim {
		case "end":
			i1, i2 = max-1, max
		case ":":
			i1, i2 = 0, max
		default:
			i1, i2 = parseRange(dim, max)
		}
	case int:
		i1, i2 = dim, dim+1
	}
	return
}

func parseRange(dim string, max int) (i1, i2 int) {
	var (
		splits = strings.Split(dim, ":")
		err    error
	)
	if i1, err = strconv.Atoi(splits[0]); err != nil {
		i1 = 0
	}
	if len(splits) == 1 {
		i2 = i1 + 1
		return
	}
	if i2, err = strconv.Atoi(splits[1]); err != nil {
		i2 = max
	}
	if i2 == i1 {
		i2 = i1 + 1
	}
	return
}

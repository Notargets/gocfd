package utils

import (
	"fmt"
	"strconv"
	"strings"
)

type R1 struct {
	Max  int
	Data []int
}

func NewR1(imax int) R1 {
	return R1{imax, nil}
}

func (r R1) Range(dimI interface{}) Index {
	var (
		i1, i2 = ParseDim(dimI, r.Max)
	)
	return NewRange(i1, i2)
}

type R2 struct {
	Ir, Jr R1
	Data   []int
}

func NewR2(imax, jmax int) R2 {
	return R2{
		NewR1(imax),
		NewR1(jmax),
		nil,
	}
}

func (r R2) Dims() (ni, nj int) {
	ni, nj = r.Ir.Max, r.Jr.Max
	return
}

func (r *R2) allocateData() {
	var (
		ni, nj = r.Dims()
	)
	if len(r.Data) == 0 {
		r.Data = make([]int, ni*nj)
	}
}

func (r *R2) Get(dimI, dimJ interface{}) (V Index) {
	var (
		rI = r.Range(dimI, dimJ)
	)
	V = make(Index, len(rI))
	r.allocateData()
	for i, ind := range rI {
		V[i] = r.Data[ind]
	}
	return
}

func (r *R2) Index() (V Index) {
	return r.Data
}

func (r *R2) Assign(dimI, dimJ interface{}, V Index) (err error) {
	r.allocateData()
	rI := r.Range(dimI, dimJ)
	if len(rI) != len(V) {
		err = fmt.Errorf("range mismatch: leftLen, rightLen = %v, %v", len(rI), len(V))
		return
	}
	for i, ind := range rI {
		r.Data[ind] = V[i]
	}
	return
}

func (r R2) Range(dimI, dimJ interface{}) (I Index) {
	var (
		i1, i2 = ParseDim(dimI, r.Ir.Max)
		j1, j2 = ParseDim(dimJ, r.Jr.Max)
		//_, nj  = r.Ir.Max, r.Jr.Max
		//ni, _  = r.Ir.Max, r.Jr.Max
		ni, nj = r.Ir.Max, r.Jr.Max
	)
	switch {
	case i1 > ni || i2 > ni:
		msg := fmt.Errorf("first dimension exceeds max, i1, i2, ni = %d, %d, %d\n", i1, i2, ni)
		panic(msg)
	case j1 > nj || j2 > nj:
		msg := fmt.Errorf("first dimension exceeds max, j1, j2, nj = %d, %d, %d\n", j1, j2, nj)
		panic(msg)
	}
	size := (i2 - i1) * (j2 - j1)
	I = NewIndex(size)
	var ind int
	for i := i1; i < i2; i++ {
		for j := j1; j < j2; j++ {
			I[ind] = j + nj*i // Column Major
			//I[ind] = i + ni*j // Row Major
			ind++
		}
	}
	return
}

type R3 struct {
	Ir, Jr, Kr R1
	Data       []int
}

func NewR3(imax, jmax, kmax int) R3 {
	return R3{
		NewR1(imax),
		NewR1(jmax),
		NewR1(kmax),
		nil,
	}
}

func (r R3) Dims() (ni, nj, nk int) {
	ni, nj, nk = r.Ir.Max, r.Jr.Max, r.Kr.Max
	return
}

func (r *R3) allocateData() {
	var (
		ni, nj, nk = r.Dims()
	)
	if len(r.Data) == 0 {
		r.Data = make([]int, ni*nj*nk)
	}
}

func (r *R3) Assign(dimI, dimJ, dimK interface{}, V Index) (err error) {
	r.allocateData()
	rI := r.Range(dimI, dimJ, dimK)
	if len(rI) != len(V) {
		err = fmt.Errorf("range mismatch: leftLen, rightLen = %v, %v", len(rI), len(V))
		return
	}
	for i, ind := range rI {
		r.Data[ind] = V[i]
	}
	return
}

func (r *R3) Index() (V Index) {
	return r.Data
}

func (r *R3) Get(dimI, dimJ, dimK interface{}) (V Index) {
	var (
		rI = r.Range(dimI, dimJ, dimK)
	)
	V = make(Index, len(rI))
	r.allocateData()
	for i, ind := range rI {
		V[i] = r.Data[ind]
	}
	return
}

func (r R3) Range(dimI, dimJ, dimK interface{}) (I Index) {
	var (
		i1, i2    = ParseDim(dimI, r.Ir.Max)
		j1, j2    = ParseDim(dimJ, r.Jr.Max)
		k1, k2    = ParseDim(dimK, r.Kr.Max)
		_, nj, nk = r.Ir.Max, r.Jr.Max, r.Kr.Max
		//ni, nj, _ = r.Ir.Max, r.Jr.Max, r.Kr.Max
	)
	size := (i2 - i1) * (j2 - j1) * (k2 - k1)
	I = NewIndex(size)
	var ind int
	for k := k1; k < k2; k++ {
		for j := j1; j < j2; j++ {
			for i := i1; i < i2; i++ {
				I[ind] = k + nk*(j+nj*i) // Column Major
				//I[ind] = i + ni*(j+nj*k) // Row Major
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
		dimC := dimInt(dim, max)
		i1, i2 = dimC, dimC+1
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
	i2 = dimInt(i2, max)
	if len(splits) == 1 {
		i2 = i1 + 1
		return
	}
	if i2, err = strconv.Atoi(splits[1]); err != nil {
		i2 = max
	}
	i2 = dimInt(i2, max)
	if i2 == i1 {
		i2 = i1 + 1
	}
	return
}

func dimInt(input, max int) (res int) {
	switch {
	case input < 0:
		// -1 means max (loop index, so max-1)
		res = max + input
	default:
		res = input
	}
	return
}

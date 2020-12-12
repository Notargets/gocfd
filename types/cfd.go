package types

import (
	"fmt"
	"sort"
	"strings"
)

//go:generate stringer -type=BCFLAG

type BCFLAG uint8

const (
	BC_None BCFLAG = iota
	BC_In
	BC_Dirichlet
	BC_Slip
	BC_Far
	BC_Wall
	BC_Cyl
	BC_Neuman
	BC_Out
	BC_IVortex
	BC_Periodic
	BC_PeriodicReversed
)

var BCNameMap = map[string]BCFLAG{
	"inflow":           BC_In,
	"in":               BC_In,
	"out":              BC_Out,
	"outflow":          BC_Out,
	"wall":             BC_Wall,
	"far":              BC_Far,
	"cyl":              BC_Cyl,
	"dirichlet":        BC_Dirichlet,
	"neuman":           BC_Neuman,
	"slip":             BC_Slip,
	"periodic":         BC_Periodic,
	"periodicreversed": BC_PeriodicReversed,
}

type BCMAP map[BCTAG][]EdgeInt // Map of BCs, key is BC tag name, e.g. "Periodic-1" or "Wall" or "Wall-top"

type BCTAG string // Tag used to name a BC consisting of a primary name ("Wall") with an optional label ("Wall-top")

func NewBCTAG(label string) (bf BCTAG) {
	bf = BCTAG(strings.ToLower(strings.Trim(label, " ")))
	return
}

func (bt BCTAG) GetFLAG() (bf BCFLAG) {
	var (
		base = string(bt)
		ind  = strings.Index(base, "-")
		ok   bool
		err  error
	)
	if ind > 0 {
		base = base[0:ind]
	}
	if bf, ok = BCNameMap[base]; !ok {
		err = fmt.Errorf("unable to find BC with base name: [%s], full tag: [%s]\n", base, string(bt))
		panic(err)
	}
	return
}

func (bt BCTAG) GetLabel() (label string) {
	var (
		base = string(bt)
		ind  = strings.Index(base, "-")
	)
	if ind > 0 {
		label = base[ind+1:]
	}
	return
}

func (bm BCMAP) AddEdges(key BCTAG, edges []EdgeInt) {
	var (
		prevInd int
		nEdges  = len(edges)
	)
	if _, ok := bm[key]; !ok {
		bm[key] = make([]EdgeInt, nEdges)
	} else {
		// Add Nedges more storage to the map
		// This will end up appending duplicate tagged BCs to a common slice
		// For instance: Periodic BCs should come in pairs, so there should be Nedges x 2 edges
		prevInd = len(bm[key])
		bsI := GrowSlice(bm[key], prevInd+nEdges)
		bm[key] = bsI.([]EdgeInt)
	}
	for i, e := range edges {
		bm[key][i+prevInd] = e
	}
}

func (bm BCMAP) Print() {
	var (
		keys = make(sort.StringSlice, len(bm))
	)
	var i int
	for key := range bm {
		keys[i] = string(key)
		i++
	}
	keys.Sort()
	for _, key := range keys {
		edges := bm[BCTAG(key)]
		for i, e := range edges {
			fmt.Printf("bc[%s][%d]=%v\n", key, i, e.GetVertices())
		}
	}
}

package Euler2D

import (
	"fmt"
	"sort"

	"github.com/notargets/gocfd/utils"
)

type VertexToElement [][2]int32 // Vertex id is the first int32, element ID is the next

func (ve VertexToElement) Len() int           { return len(ve) }
func (ve VertexToElement) Swap(i, j int)      { ve[i], ve[j] = ve[j], ve[i] }
func (ve VertexToElement) Less(i, j int) bool { return ve[i][0] < ve[j][0] }
func (ve VertexToElement) Sort()              { sort.Sort(ve) }

func NewVertexToElement(EtoV utils.Matrix) (VtoE VertexToElement) {
	var (
		Kmax, Nverts = EtoV.Dims()
	)
	if Nverts != 3 {
		msg := fmt.Errorf("EtoV should have dimensions [Kmax,3] was [%d,%d]", Kmax, Nverts)
		panic(msg)
	}
	VtoE = make(VertexToElement, Kmax*3)
	var ii int
	for k := 0; k < Kmax; k++ {
		for i := 0; i < 3; i++ {
			VtoE[ii] = [2]int32{int32(EtoV.At(k, i)), int32(k)}
			ii++
		}
	}
	VtoE.Sort()
	return
}

type ScalarDissipation struct {
	EpsVertex    [][]float64    // Sharded Nverts, Aggregated (Max) of epsilon surrounding each vertex
	Epsilon      []utils.Matrix // Sharded Np x Kmax, Interpolated from element vertices
	DissX, DissY []utils.Matrix // Sharded Np x Kmax, Dissipation Flux
	RDiss        []utils.Matrix // Sharded Np x Kmax, Dissipation Added to Residual
}

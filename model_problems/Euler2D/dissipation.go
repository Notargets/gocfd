package Euler2D

import (
	"fmt"
	"sort"

	"github.com/notargets/gocfd/DG2D"

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

func (ve VertexToElement) Shard(NPar int) (veSharded []VertexToElement) {
	var (
		lve              = len(ve)
		VertexPartitions = NewPartitionMap(NPar, lve) // This has to be re-done to honor vertex grouping
		ib               int
		vNum             int32
	)
	veSharded = make([]VertexToElement, NPar)
	approxBucketSize := VertexPartitions.GetBucketDimension(0)
	for np := 0; np < NPar; np++ {
		for i := 0; i < approxBucketSize; i++ {
			//veSharded[np][i] = ve[ib]
			veSharded[np] = append(veSharded[np], ve[ib])
			ib++
			//fmt.Printf("i,ib,np = %d,%d,%d\n", i, ib, np)
			if ib == lve {
				return
			}
		}
		vNum = ve[ib][0]
		for ib < lve && ve[ib][0] == vNum {
			//fmt.Printf("ib,np,vNum,ve[ib][0] = %d,%d,%d,%d\n", ib, np, vNum, ve[ib][0])
			// Continue the shard until the group is done
			veSharded[np] = append(veSharded[np], ve[ib])
			ib++
			if ib == lve {
				return
			}
		}
	}
	return
}

type ScalarDissipation struct {
	VtoE             []VertexToElement // Sharded vertex to element map, [2] is [vertID, ElementID]
	Epsilon          []utils.Matrix    // Sharded Np x Kmax, Interpolated from element vertices
	DissX, DissY     [][4]utils.Matrix // Sharded Np x Kmax, Dissipation Flux
	RDiss            [][4]utils.Matrix // Sharded Np x Kmax, Dissipation Added to Residual
	EpsVertex        [][]float64       // Sharded Nverts, Aggregated (Max) of epsilon surrounding each vertex
	VertexPartitions *PartitionMap
}

func NewScalarDissipation(NVerts int, EToV utils.Matrix, pm *PartitionMap, el *DG2D.LagrangeElement2D) (sd *ScalarDissipation) {
	var (
		NPar = pm.ParallelDegree
		Np   = el.Np
	)
	sd = &ScalarDissipation{
		EpsVertex: make([][]float64, NPar),
		Epsilon:   make([]utils.Matrix, NPar),
		DissX:     make([][4]utils.Matrix, NPar),
		DissY:     make([][4]utils.Matrix, NPar),
		RDiss:     make([][4]utils.Matrix, NPar),
		VtoE:      NewVertexToElement(EToV).Shard(NPar),
	}
	for np := 0; np < NPar; np++ {
		sd.EpsVertex[np] = make([]float64, NVerts)
		Kmax := pm.GetBucketDimension(np)
		sd.Epsilon[np] = utils.NewMatrix(Np, Kmax)
		for n := 0; n < 4; n++ {
			sd.DissX[np][n] = utils.NewMatrix(Np, Kmax)
			sd.DissY[np][n] = utils.NewMatrix(Np, Kmax)
			sd.RDiss[np][n] = utils.NewMatrix(Np, Kmax)
		}
	}
	return
}

package tetelement

import (
	"fmt"
	"github.com/notargets/gocfd/DG3D/mesh"
	"github.com/notargets/gocfd/DG3D/mesh/readers"
	"github.com/notargets/gocfd/DG3D/tetrahedra/gonudg"
)

type Element3D struct {
	K    int   // Number of elements
	EToP []int // Element to Parallel partitions
	*gonudg.DG3D
	Mesh           *mesh.Mesh
	SplitElement3D []*Element3D
}

// NewElement3D creates an Element3D from a mesh file
func NewElement3D(order int, meshFile string) (el *Element3D, err error) {
	// Read mesh file
	m, err := readers.ReadMeshFile(meshFile)
	if err != nil {
		return nil, err
	}

	// Use the common constructor
	return NewElement3DFromMesh(order, m)
}

// NewElement3DFromMesh creates an Element3D from an existing mesh
func NewElement3DFromMesh(order int, m *mesh.Mesh) (el *Element3D, err error) {
	el = &Element3D{
		Mesh: m,
	}
	nverts := len(m.Vertices)
	VX, VY, VZ := make([]float64, nverts), make([]float64, nverts), make([]float64, nverts)
	el.DG3D, err = gonudg.NewDG3D(order, VX, VY, VZ, m.EtoV)
	if err != nil {
		fmt.Println("Unable to convert mesh to DG3D: ", err)
		panic(err)
	}

	// Split mesh by partition if EToP is present
	if el.EToP != nil {
		if err = el.SplitByPartition(); err != nil {
			return nil, fmt.Errorf("failed to split mesh by partition: %v", err)
		}
	}

	return el, nil
}

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
	BCMaps         *BCFaceMap // Boundary condition mapping
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
		K:    m.NumElements,
	}
	nverts := len(m.Vertices)
	VX, VY, VZ := make([]float64, nverts), make([]float64, nverts), make([]float64, nverts)

	// Copy vertex coordinates
	for i := 0; i < nverts; i++ {
		VX[i] = m.Vertices[i][0]
		VY[i] = m.Vertices[i][1]
		VZ[i] = m.Vertices[i][2]
	}

	el.DG3D, err = gonudg.NewDG3D(order, VX, VY, VZ, m.EtoV)
	if err != nil {
		fmt.Println("Unable to convert mesh to DG3D: ", err)
		panic(err)
	}

	// Set up connectivity on DG3D object before BuildMaps3D
	el.DG3D.EToE = m.EToE
	el.DG3D.EToF = m.EToF

	// Build DG3D maps (needs connectivity to identify boundary faces)
	el.DG3D.BuildMaps3D()

	// Build boundary condition maps
	if err = el.BuildBCMaps(); err != nil {
		return nil, fmt.Errorf("failed to build BC maps: %v", err)
	}

	// Copy EToP if present
	el.EToP = m.EToP

	// Split mesh by partition if EToP is present
	if el.EToP != nil {
		if err = el.SplitByPartition(); err != nil {
			return nil, fmt.Errorf("failed to split mesh by partition: %v", err)
		}
	}

	return el, nil
}

// GetBoundaryFaces returns the element and face indices for all boundary faces
func (el *Element3D) GetBoundaryFaces() (elements []int, faces []int) {
	elements = make([]int, 0)
	faces = make([]int, 0)

	if el.EToE == nil {
		return
	}

	for k := 0; k < el.K; k++ {
		for f := 0; f < el.Nfaces; f++ {
			// Check if this is a boundary face
			if el.EToE[k][f] == k && el.EToF[k][f] == f {
				elements = append(elements, k)
				faces = append(faces, f)
			}
		}
	}

	return
}

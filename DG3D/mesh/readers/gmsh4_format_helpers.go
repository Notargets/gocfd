package readers

import (
	"fmt"
	"github.com/notargets/gocfd/DG3D/mesh"
	"strings"
)

// Gmsh4TestBuilder helps build Gmsh 4.1 format test files
type Gmsh4TestBuilder struct {
	tm *mesh.TestMeshes
}

// NewGmsh4TestBuilder creates a new builder with standard test meshes
func NewGmsh4TestBuilder() *Gmsh4TestBuilder {
	return &Gmsh4TestBuilder{
		tm: mesh.GetStandardTestMeshes(),
	}
}

// BuildMixedElementTest creates a Gmsh 4.1 file with mixed element types
func (b *Gmsh4TestBuilder) BuildMixedElementTest() string {
	mesh := b.tm.MixedMesh
	return b.BuildFromCompleteMesh(&mesh)
}

// BuildTwoTetTest creates a Gmsh 4.1 file with two tetrahedra
func (b *Gmsh4TestBuilder) BuildTwoTetTest() string {
	mesh := b.tm.TwoTetMesh
	return b.BuildFromCompleteMesh(&mesh)
}

// BuildFromCompleteMesh creates a complete Gmsh 4.1 format file from a CompleteMesh
func (b *Gmsh4TestBuilder) BuildFromCompleteMesh(mesh *mesh.CompleteMesh) string {
	var sections []string

	// Header
	sections = append(sections, b.buildHeader())

	// Entities
	sections = append(sections, b.buildEntities(mesh))

	// Nodes
	sections = append(sections, b.buildNodes(mesh))

	// Elements
	sections = append(sections, b.buildElements(mesh))

	return strings.Join(sections, "\n")
}

func (b *Gmsh4TestBuilder) buildHeader() string {
	return `$MeshFormat
4.1 0 8
$EndMeshFormat`
}

func (b *Gmsh4TestBuilder) buildEntities(mesh *mesh.CompleteMesh) string {
	// For simplicity, create one volume entity
	minX, minY, minZ := mesh.BoundingBox[0][0], mesh.BoundingBox[0][1], mesh.BoundingBox[0][2]
	maxX, maxY, maxZ := mesh.BoundingBox[1][0], mesh.BoundingBox[1][1], mesh.BoundingBox[1][2]

	return fmt.Sprintf(`$Entities
0 0 0 1
1 %f %f %f %f %f %f 0
$EndEntities`, minX, minY, minZ, maxX, maxY, maxZ)
}

func (b *Gmsh4TestBuilder) buildNodes(mesh *mesh.CompleteMesh) string {
	numNodes := len(mesh.Nodes.Nodes)

	var lines []string
	lines = append(lines, "$Nodes")
	lines = append(lines, fmt.Sprintf("1 %d 1 %d", numNodes, numNodes))
	lines = append(lines, fmt.Sprintf("%d 1 0 %d", mesh.Dimension, numNodes))

	// Node tags
	for i := 1; i <= numNodes; i++ {
		lines = append(lines, fmt.Sprintf("%d", i))
	}

	// Node coordinates
	for i := 0; i < numNodes; i++ {
		coords := mesh.Nodes.Nodes[i]
		lines = append(lines, fmt.Sprintf("%f %f %f", coords[0], coords[1], coords[2]))
	}

	lines = append(lines, "$EndNodes")
	return strings.Join(lines, "\n")
}

func (b *Gmsh4TestBuilder) buildElements(mesh *mesh.CompleteMesh) string {
	// Count total elements and blocks
	totalElements := 0
	numBlocks := 0
	for _, elemSet := range mesh.Elements {
		if len(elemSet.Elements) > 0 {
			numBlocks++
			totalElements += len(elemSet.Elements)
		}
	}

	var lines []string
	lines = append(lines, "$Elements")
	lines = append(lines, fmt.Sprintf("%d %d 1 %d", numBlocks, totalElements, totalElements))

	elemTag := 1
	for _, elemSet := range mesh.Elements {
		if len(elemSet.Elements) == 0 {
			continue
		}

		// Block header: entityDim entityTag elementType numElements
		lines = append(lines, fmt.Sprintf("3 1 %d %d",
			elementTypeToGmsh4[elemSet.Type], len(elemSet.Elements)))

		// Elements
		for _, elem := range elemSet.Elements {
			nodeIDs := make([]string, len(elem))
			for i, nodeName := range elem {
				nodeIDs[i] = fmt.Sprintf("%d", mesh.Nodes.NodeIDMap[nodeName])
			}
			lines = append(lines, fmt.Sprintf("%d %s", elemTag, strings.Join(nodeIDs, " ")))
			elemTag++
		}
	}

	lines = append(lines, "$EndElements")
	return strings.Join(lines, "\n")
}

// Helper to convert our ElementType to Gmsh element type number
var elementTypeToGmsh4 = map[mesh.ElementType]int{
	mesh.Point:      15,
	mesh.Line:       1,
	mesh.Line3:      8,
	mesh.Triangle:   2,
	mesh.Triangle6:  9,
	mesh.Triangle9:  20,
	mesh.Triangle10: 21,
	mesh.Quad:       3,
	mesh.Quad8:      16,
	mesh.Quad9:      10,
	mesh.Tet:        4,
	mesh.Tet10:      11,
	mesh.Hex:        5,
	mesh.Hex20:      17,
	mesh.Hex27:      12,
	mesh.Prism:      6,
	mesh.Prism15:    18,
	mesh.Prism18:    13,
	mesh.Pyramid:    7,
	mesh.Pyramid13:  19,
	mesh.Pyramid14:  14,
}

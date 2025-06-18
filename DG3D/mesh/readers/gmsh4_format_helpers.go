package readers

import (
	"fmt"
	"github.com/notargets/gocfd/utils"
	"strings"
)

// Gmsh4TestBuilder helps build Gmsh 4.1 format test files
type Gmsh4TestBuilder struct {
	tm *utils.TestMeshes
}

// NewGmsh4TestBuilder creates a new builder with standard test meshes
func NewGmsh4TestBuilder() *Gmsh4TestBuilder {
	return &Gmsh4TestBuilder{
		tm: utils.GetStandardTestMeshes(),
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
func (b *Gmsh4TestBuilder) BuildFromCompleteMesh(mesh *utils.CompleteMesh) string {
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

func (b *Gmsh4TestBuilder) buildEntities(mesh *utils.CompleteMesh) string {
	// For simplicity, create one volume entity
	minX, minY, minZ := mesh.BoundingBox[0][0], mesh.BoundingBox[0][1], mesh.BoundingBox[0][2]
	maxX, maxY, maxZ := mesh.BoundingBox[1][0], mesh.BoundingBox[1][1], mesh.BoundingBox[1][2]

	return fmt.Sprintf(`$Entities
0 0 0 1
1 %f %f %f %f %f %f 0
$EndEntities`, minX, minY, minZ, maxX, maxY, maxZ)
}

func (b *Gmsh4TestBuilder) buildNodes(mesh *utils.CompleteMesh) string {
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

func (b *Gmsh4TestBuilder) buildElements(mesh *utils.CompleteMesh) string {
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
var elementTypeToGmsh4 = map[utils.ElementType]int{
	utils.Point:      15,
	utils.Line:       1,
	utils.Line3:      8,
	utils.Triangle:   2,
	utils.Triangle6:  9,
	utils.Triangle9:  20,
	utils.Triangle10: 21,
	utils.Quad:       3,
	utils.Quad8:      16,
	utils.Quad9:      10,
	utils.Tet:        4,
	utils.Tet10:      11,
	utils.Hex:        5,
	utils.Hex20:      17,
	utils.Hex27:      12,
	utils.Prism:      6,
	utils.Prism15:    18,
	utils.Prism18:    13,
	utils.Pyramid:    7,
	utils.Pyramid13:  19,
	utils.Pyramid14:  14,
}

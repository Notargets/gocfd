package readers

import (
	"fmt"
	"github.com/notargets/gocfd/utils"
	"strings"
)

// Gmsh22TestBuilder helps build Gmsh 2.2 format test files
type Gmsh22TestBuilder struct {
	tm *utils.TestMeshes
}

// NewGmsh22TestBuilder creates a new builder with standard test meshes
func NewGmsh22TestBuilder() *Gmsh22TestBuilder {
	return &Gmsh22TestBuilder{
		tm: utils.GetStandardTestMeshes(),
	}
}

// BuildMixedElementTest creates a Gmsh 2.2 file with mixed element types
func (b *Gmsh22TestBuilder) BuildMixedElementTest() string {
	mesh := b.tm.MixedMesh
	return b.BuildFromCompleteMesh(&mesh)
}

// BuildTwoTetTest creates a Gmsh 2.2 file with two tetrahedra
func (b *Gmsh22TestBuilder) BuildTwoTetTest() string {
	mesh := b.tm.TwoTetMesh
	return b.BuildFromCompleteMesh(&mesh)
}

// BuildCubeTest creates a Gmsh 2.2 file with the cube mesh
func (b *Gmsh22TestBuilder) BuildCubeTest() string {
	mesh := b.tm.CubeMesh
	return b.BuildFromCompleteMesh(&mesh)
}

// BuildFromCompleteMesh creates a complete Gmsh 2.2 format file from a CompleteMesh
func (b *Gmsh22TestBuilder) BuildFromCompleteMesh(mesh *utils.CompleteMesh) string {
	var sections []string

	// Header
	sections = append(sections, b.buildHeader())

	// Nodes
	sections = append(sections, b.buildNodes(mesh))

	// Elements
	sections = append(sections, b.buildElements(mesh))

	return strings.Join(sections, "\n")
}

func (b *Gmsh22TestBuilder) buildHeader() string {
	return `$MeshFormat
2.2 0 8
$EndMeshFormat`
}

func (b *Gmsh22TestBuilder) buildNodes(mesh *utils.CompleteMesh) string {
	numNodes := len(mesh.Nodes.Nodes)

	var lines []string
	lines = append(lines, "$Nodes")
	lines = append(lines, fmt.Sprintf("%d", numNodes))

	// Node lines: id x y z
	for i := 0; i < numNodes; i++ {
		nodeID := i + 1 // 1-based
		coords := mesh.Nodes.Nodes[i]
		lines = append(lines, fmt.Sprintf("%d %f %f %f", nodeID, coords[0], coords[1], coords[2]))
	}

	lines = append(lines, "$EndNodes")
	return strings.Join(lines, "\n")
}

func (b *Gmsh22TestBuilder) buildElements(msh *utils.CompleteMesh) string {
	// Count total elements
	totalElements := 0
	for _, elemSet := range msh.Elements {
		totalElements += len(elemSet.Elements)
	}

	var lines []string
	lines = append(lines, "$Elements")
	lines = append(lines, fmt.Sprintf("%d", totalElements))

	elemID := 1
	for _, elemSet := range msh.Elements {
		gmshType := elementTypeToGmsh22[elemSet.Type]

		for i, elem := range elemSet.Elements {
			// Get properties
			props := utils.ElementProps{}
			if i < len(elemSet.Properties) {
				props = elemSet.Properties[i]
			}

			// Build tags
			var tags []int
			if props.PhysicalTag > 0 {
				tags = append(tags, props.PhysicalTag)
			}
			if props.GeometricTag > 0 {
				tags = append(tags, props.GeometricTag)
			}

			// Default to minimal tags if none specified
			if len(tags) == 0 {
				tags = []int{0}
			}

			// Convert node names to IDs
			nodeIDs := make([]string, len(elem))
			for j, nodeName := range elem {
				nodeIDs[j] = fmt.Sprintf("%d", msh.Nodes.NodeIDMap[nodeName])
			}

			// Format: elem-id elem-type num-tags tag1 tag2 ... node1 node2 ...
			line := fmt.Sprintf("%d %d %d", elemID, gmshType, len(tags))
			for _, tag := range tags {
				line += fmt.Sprintf(" %d", tag)
			}
			line += " " + strings.Join(nodeIDs, " ")

			lines = append(lines, line)
			elemID++
		}
	}

	lines = append(lines, "$EndElements")
	return strings.Join(lines, "\n")
}

// Helper to convert our ElementType to Gmsh 2.2 element type number
var elementTypeToGmsh22 = map[utils.ElementType]int{
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

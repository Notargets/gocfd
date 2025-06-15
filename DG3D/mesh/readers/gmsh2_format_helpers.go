package readers

import (
	"fmt"
	"github.com/notargets/gocfd/DG3D/mesh"
	"strings"
)

// Gmsh22TestBuilder helps build Gmsh 2.2 format test files
type Gmsh22TestBuilder struct {
	tm *mesh.TestMeshes
}

// NewGmsh22TestBuilder creates a new builder with standard test meshes
func NewGmsh22TestBuilder() *Gmsh22TestBuilder {
	return &Gmsh22TestBuilder{
		tm: mesh.GetStandardTestMeshes(),
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
func (b *Gmsh22TestBuilder) BuildFromCompleteMesh(mesh *mesh.CompleteMesh) string {
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

func (b *Gmsh22TestBuilder) buildNodes(mesh *mesh.CompleteMesh) string {
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

func (b *Gmsh22TestBuilder) buildElements(msh *mesh.CompleteMesh) string {
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
			props := mesh.ElementProps{}
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
var elementTypeToGmsh22 = map[mesh.ElementType]int{
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

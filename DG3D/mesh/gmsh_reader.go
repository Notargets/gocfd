package mesh

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"
)

// ReadGmsh reads a Gmsh format file (version 2.2)
func ReadGmsh(filename string) (*Mesh, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	mesh := NewMesh()
	scanner := bufio.NewScanner(file)

	var version string

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())

		switch line {
		case "$MeshFormat":
			// Read format version
			scanner.Scan()
			parts := strings.Fields(scanner.Text())
			if len(parts) > 0 {
				version = parts[0]
			}
			// Skip until $EndMeshFormat
			for scanner.Scan() {
				if strings.TrimSpace(scanner.Text()) == "$EndMeshFormat" {
					break
				}
			}

		case "$Nodes":
			scanner.Scan()
			numNodes, _ := strconv.Atoi(scanner.Text())
			mesh.Vertices = make([][]float64, numNodes)

			for i := 0; i < numNodes; i++ {
				scanner.Scan()
				fields := strings.Fields(scanner.Text())
				if len(fields) >= 4 {
					id, _ := strconv.Atoi(fields[0])
					x, _ := strconv.ParseFloat(fields[1], 64)
					y, _ := strconv.ParseFloat(fields[2], 64)
					z, _ := strconv.ParseFloat(fields[3], 64)
					if id >= 1 && id <= numNodes {
						mesh.Vertices[id-1] = []float64{x, y, z}
					}
				}
			}
			// Skip $EndNodes
			scanner.Scan()

		case "$Elements":
			scanner.Scan()
			numElems, _ := strconv.Atoi(scanner.Text())

			// First pass: count 3D elements
			elemData := [][]string{}
			for i := 0; i < numElems; i++ {
				scanner.Scan()
				elemData = append(elemData, strings.Fields(scanner.Text()))
			}

			// Second pass: extract 3D elements only
			for _, fields := range elemData {
				if len(fields) < 3 {
					continue
				}

				elemType, _ := strconv.Atoi(fields[1])

				// Only process 3D elements
				var etype ElementType
				var numNodes int
				validElem := true

				switch elemType {
				case 4: // 4-node tet
					etype = Tet
					numNodes = 4
				case 5: // 8-node hex
					etype = Hex
					numNodes = 8
				case 6: // 6-node prism
					etype = Prism
					numNodes = 6
				case 7: // 5-node pyramid
					etype = Pyramid
					numNodes = 5
				case 11: // 10-node tet (2nd order)
					etype = Tet
					numNodes = 4 // Only use corner nodes
				case 12: // 27-node hex (2nd order)
					etype = Hex
					numNodes = 8 // Only use corner nodes
				case 13: // 18-node prism (2nd order)
					etype = Prism
					numNodes = 6 // Only use corner nodes
				case 14: // 14-node pyramid (2nd order)
					etype = Pyramid
					numNodes = 5 // Only use corner nodes
				default:
					validElem = false
				}

				if validElem {
					numTags, _ := strconv.Atoi(fields[2])
					physTag := 0
					if numTags > 0 && len(fields) > 3 {
						physTag, _ = strconv.Atoi(fields[3])
					}

					// Check if we have enough fields
					offset := 3 + numTags
					if len(fields) >= offset+numNodes {
						// Read vertices (convert to 0-indexed)
						verts := make([]int, numNodes)
						for j := 0; j < numNodes; j++ {
							v, _ := strconv.Atoi(fields[offset+j])
							verts[j] = v - 1
						}

						mesh.Elements = append(mesh.Elements, verts)
						mesh.ElementTypes = append(mesh.ElementTypes, etype)
						mesh.ElementTags = append(mesh.ElementTags, physTag)
					}
				}
			}
			// Skip $EndElements
			scanner.Scan()

		case "$PhysicalNames":
			scanner.Scan()
			numPhysical, _ := strconv.Atoi(scanner.Text())

			for i := 0; i < numPhysical; i++ {
				scanner.Scan()
				fields := strings.Fields(scanner.Text())
				if len(fields) >= 3 {
					// dim, _ := strconv.Atoi(fields[0])
					tag, _ := strconv.Atoi(fields[1])
					name := strings.Trim(fields[2], "\"")
					mesh.BoundaryTags[tag] = name
				}
			}
			// Skip $EndPhysicalNames
			scanner.Scan()
		}
	}

	// Handle version 4.1 format if needed
	if strings.HasPrefix(version, "4") {
		return nil, fmt.Errorf("Gmsh format version %s not yet implemented", version)
	}

	mesh.NumElements = len(mesh.Elements)
	mesh.NumVertices = len(mesh.Vertices)
	mesh.BuildConnectivity()

	return mesh, nil
}

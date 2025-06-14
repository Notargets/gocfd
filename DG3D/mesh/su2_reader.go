package mesh

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"
)

// ReadSU2 reads an SU2 native format file
func ReadSU2(filename string) (*Mesh, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	mesh := NewMesh()
	scanner := bufio.NewScanner(file)

	var ndime int

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())

		// Skip comments
		if strings.HasPrefix(line, "%") || line == "" {
			continue
		}

		if strings.HasPrefix(line, "NDIME=") {
			fmt.Sscanf(line, "NDIME=%d", &ndime)
			if ndime != 3 {
				return nil, fmt.Errorf("only 3D meshes are supported, got NDIME=%d", ndime)
			}

		} else if strings.HasPrefix(line, "NELEM=") {
			var nelem int
			fmt.Sscanf(line, "NELEM=%d", &nelem)

			// Pre-allocate slices
			mesh.Elements = make([][]int, 0, nelem)
			mesh.ElementTypes = make([]ElementType, 0, nelem)
			mesh.ElementTags = make([]int, 0, nelem)

			// Read elements
			for i := 0; i < nelem; i++ {
				scanner.Scan()
				fields := strings.Fields(scanner.Text())

				if len(fields) < 2 {
					continue
				}

				su2Type, _ := strconv.Atoi(fields[0])

				// Map SU2 element types to our types
				var etype ElementType
				var numNodes int
				validElem := true

				switch su2Type {
				case 3: // Line
					validElem = false // Skip 1D elements
				case 5: // Triangle
					validElem = false // Skip 2D elements
				case 9: // Quad
					validElem = false // Skip 2D elements
				case 10: // Tet
					etype = Tet
					numNodes = 4
				case 12: // Hex
					etype = Hex
					numNodes = 8
				case 13: // Prism
					etype = Prism
					numNodes = 6
				case 14: // Pyramid
					etype = Pyramid
					numNodes = 5
				default:
					validElem = false
				}

				if validElem && len(fields) >= numNodes+2 {
					// Read vertices
					verts := make([]int, numNodes)
					for j := 0; j < numNodes; j++ {
						verts[j], _ = strconv.Atoi(fields[1+j])
					}

					// Element ID is last field
					// elemID, _ := strconv.Atoi(fields[len(fields)-1])

					mesh.Elements = append(mesh.Elements, verts)
					mesh.ElementTypes = append(mesh.ElementTypes, etype)
					mesh.ElementTags = append(mesh.ElementTags, 0) // Default tag
				}
			}

		} else if strings.HasPrefix(line, "NPOIN=") {
			var npoin int
			fmt.Sscanf(line, "NPOIN=%d", &npoin)

			mesh.Vertices = make([][]float64, npoin)

			for i := 0; i < npoin; i++ {
				scanner.Scan()
				fields := strings.Fields(scanner.Text())

				if len(fields) >= ndime+1 {
					coords := make([]float64, 3)
					for j := 0; j < ndime; j++ {
						coords[j], _ = strconv.ParseFloat(fields[j], 64)
					}

					// Point ID is last field
					ptID, _ := strconv.Atoi(fields[len(fields)-1])
					if ptID >= 0 && ptID < npoin {
						mesh.Vertices[ptID] = coords
					}
				}
			}

		} else if strings.HasPrefix(line, "NMARK=") {
			var nmark int
			fmt.Sscanf(line, "NMARK=%d", &nmark)

			// Read boundary markers
			for i := 0; i < nmark; i++ {
				scanner.Scan()
				markerLine := strings.TrimSpace(scanner.Text())

				if strings.HasPrefix(markerLine, "MARKER_TAG=") {
					tagName := strings.TrimPrefix(markerLine, "MARKER_TAG=")
					tagName = strings.Trim(tagName, " ")

					// Read number of marker elements
					scanner.Scan()
					elemLine := strings.TrimSpace(scanner.Text())
					var nMarkerElems int
					fmt.Sscanf(elemLine, "MARKER_ELEMS=%d", &nMarkerElems)

					// Store boundary tag
					mesh.BoundaryTags[i] = tagName

					// Skip marker elements for now (could be extended to read them)
					for j := 0; j < nMarkerElems; j++ {
						scanner.Scan()
					}
				}
			}
		}
	}

	mesh.NumElements = len(mesh.Elements)
	mesh.NumVertices = len(mesh.Vertices)
	mesh.BuildConnectivity()

	return mesh, nil
}

// getNumNodesSU2 returns the number of nodes for an SU2 element type
func getNumNodesSU2(su2Type int) int {
	switch su2Type {
	case 3:
		return 2 // Line
	case 5:
		return 3 // Triangle
	case 9:
		return 4 // Quad
	case 10:
		return 4 // Tet
	case 12:
		return 8 // Hex
	case 13:
		return 6 // Prism
	case 14:
		return 5 // Pyramid
	default:
		return 0
	}
}

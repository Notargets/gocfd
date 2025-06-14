package mesh

import (
	"bufio"
	"os"
	"strconv"
	"strings"
)

// ReadGambitNeutral reads a Gambit neutral file
func ReadGambitNeutral(filename string) (*Mesh, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	mesh := NewMesh()
	scanner := bufio.NewScanner(file)

	var numnp, nelem int

	// Read until we find the problem size parameters
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())

		// Look for the problem size header line
		if strings.Contains(line, "NUMNP") && strings.Contains(line, "NELEM") {
			// Next line contains the actual values
			if scanner.Scan() {
				values := strings.Fields(scanner.Text())
				if len(values) >= 2 {
					numnp, _ = strconv.Atoi(values[0])
					nelem, _ = strconv.Atoi(values[1])
				}
			}
			break
		}
	}

	// Continue reading sections
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())

		if line == "ENDOFSECTION" {
			continue
		}

		if strings.Contains(line, "NODAL COORDINATES") {
			// Read nodes
			mesh.Vertices = make([][]float64, numnp)
			nodeCount := 0

			for scanner.Scan() {
				line = strings.TrimSpace(scanner.Text())
				if line == "ENDOFSECTION" {
					break
				}

				fields := strings.Fields(line)
				if len(fields) >= 4 {
					id, _ := strconv.Atoi(fields[0])
					x, _ := strconv.ParseFloat(fields[1], 64)
					y, _ := strconv.ParseFloat(fields[2], 64)
					z, _ := strconv.ParseFloat(fields[3], 64)

					if id >= 1 && id <= numnp {
						mesh.Vertices[id-1] = []float64{x, y, z}
						nodeCount++
					}
				}
			}

		} else if strings.Contains(line, "ELEMENTS/CELLS") {
			// Read elements
			mesh.Elements = make([][]int, 0, nelem)
			mesh.ElementTypes = make([]ElementType, 0, nelem)
			mesh.ElementTags = make([][]int, 0, nelem)

			for scanner.Scan() {
				line = strings.TrimSpace(scanner.Text())
				if line == "ENDOFSECTION" {
					break
				}

				fields := strings.Fields(line)
				if len(fields) >= 3 {
					// Format: NE NTYPE NDP NODE1 NODE2 ...
					// id, _ := strconv.Atoi(fields[0]) // Element ID (not used for indexing)
					elemType, _ := strconv.Atoi(fields[1])
					numNodes, _ := strconv.Atoi(fields[2])

					// Only process 3D elements
					var etype ElementType
					validElem := true
					switch elemType {
					case 1: // Edge - skip
						validElem = false
					case 2: // Quad - skip 2D
						validElem = false
					case 3: // Triangle - skip 2D
						validElem = false
					case 4: // Brick
						etype = Hex
					case 5: // Wedge/Prism
						etype = Prism
					case 6: // Tetrahedron
						etype = Tet
					case 7: // Pyramid
						etype = Pyramid
					default:
						validElem = false
					}

					if validElem && len(fields) >= 3+numNodes {
						// Read vertices (convert to 0-indexed)
						verts := make([]int, numNodes)
						for j := 0; j < numNodes; j++ {
							v, _ := strconv.Atoi(fields[3+j])
							verts[j] = v - 1
						}

						mesh.Elements = append(mesh.Elements, verts)
						mesh.ElementTypes = append(mesh.ElementTypes, etype)
						mesh.ElementTags = append(mesh.ElementTags,
							[]int{0}) // Default tag
					}
				}
			}

		} else if strings.Contains(line, "ELEMENT GROUP") {
			// Read element group information
			// Format: GROUP: NGP ELEMENTS: NELGP MATERIAL: MTYP NFLAGS: NFLAGS
			if strings.HasPrefix(line, "GROUP:") {
				parts := strings.Split(line, " ")
				var groupID []int
				var numElems int
				for i := 0; i < len(parts)-1; i++ {
					if parts[i] == "GROUP:" && i+1 < len(parts) {
						iii, _ := strconv.Atoi(parts[i+1])
						groupID = []int{iii}
					}
					if parts[i] == "ELEMENTS:" && i+1 < len(parts) {
						numElems, _ = strconv.Atoi(parts[i+1])
					}
				}

				// Skip entity name
				scanner.Scan()

				// Skip flags
				scanner.Scan()

				// Read element IDs in this group
				elementsRead := 0
				for scanner.Scan() && elementsRead < numElems {
					line = strings.TrimSpace(scanner.Text())
					if line == "ENDOFSECTION" {
						break
					}

					fields := strings.Fields(line)
					for _, field := range fields {
						elemID, _ := strconv.Atoi(field)
						if elemID > 0 && elemID <= len(mesh.Elements) {
							// Elements are 1-indexed in file
							mesh.ElementTags[elemID-1] = groupID
						}
						elementsRead++
						if elementsRead >= numElems {
							break
						}
					}
				}
			}

		} else if strings.Contains(line, "BOUNDARY CONDITIONS") {
			// Read boundary conditions
			scanner.Scan() // Read BC name and info
			bcLine := strings.TrimSpace(scanner.Text())
			parts := strings.Fields(bcLine)
			if len(parts) >= 2 {
				bcName := parts[0]
				// Could parse ITYPE, NENTRY, etc. if needed

				// For now, just store the BC name
				// You could extend this to read the actual boundary elements/nodes
				mesh.BoundaryTags[len(mesh.BoundaryTags)] = bcName
			}
		}
	}

	mesh.NumElements = len(mesh.Elements)
	mesh.NumVertices = len(mesh.Vertices)
	mesh.BuildConnectivity()

	return mesh, nil
}

package readers

import (
	"bufio"
	"fmt"
	"github.com/notargets/gocfd/DG3D/mesh"
	"os"
	"strconv"
	"strings"
)

// ReadGambitNeutral reads a Gambit neutral file (.neu)
func ReadGambitNeutral(filename string) (*mesh.Mesh, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	msh := mesh.NewMesh()
	scanner := bufio.NewScanner(file)

	// Control variables from header
	var numnp, nelem, ngrps, nbsets int

	// Read control info section
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if strings.Contains(line, "NUMNP") && strings.Contains(line, "NELEM") {
			// Next line contains the actual values
			if !scanner.Scan() {
				return nil, fmt.Errorf("unexpected EOF after control header")
			}
			values := strings.Fields(scanner.Text())
			if len(values) >= 4 {
				numnp, _ = strconv.Atoi(values[0])  // Number of nodes
				nelem, _ = strconv.Atoi(values[1])  // Number of elements
				ngrps, _ = strconv.Atoi(values[2])  // Number of element groups
				nbsets, _ = strconv.Atoi(values[3]) // Number of boundary sets
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
			msh.Vertices = make([][]float64, numnp)
			msh.NodeIDMap = make(map[int]int)
			msh.NodeArrayMap = make(map[int]int)

			for i := 0; i < numnp; i++ {
				if !scanner.Scan() {
					return nil, fmt.Errorf("unexpected EOF reading nodes")
				}

				fields := strings.Fields(scanner.Text())
				if len(fields) >= 4 {
					nodeID, _ := strconv.Atoi(fields[0])
					x, _ := strconv.ParseFloat(fields[1], 64)
					y, _ := strconv.ParseFloat(fields[2], 64)
					z, _ := strconv.ParseFloat(fields[3], 64)

					// Gambit uses 1-based node IDs
					idx := nodeID - 1
					if idx >= 0 && idx < numnp {
						msh.Vertices[idx] = []float64{x, y, z}
						msh.NodeIDMap[nodeID] = idx
						msh.NodeArrayMap[idx] = nodeID
					}
				}
			}

		} else if strings.Contains(line, "ELEMENTS/CELLS") {
			// Pre-allocate element arrays
			msh.EtoV = make([][]int, 0, nelem)
			msh.ElementTypes = make([]mesh.ElementType, 0, nelem)
			msh.ElementTags = make([][]int, 0, nelem)

			// Read elements
			for i := 0; i < nelem; i++ {
				if !scanner.Scan() {
					return nil, fmt.Errorf("unexpected EOF reading elements")
				}

				fields := strings.Fields(scanner.Text())
				if len(fields) >= 3 {
					// elemID, _ := strconv.Atoi(fields[0])
					gambitType, _ := strconv.Atoi(fields[1])
					numNodes, _ := strconv.Atoi(fields[2])

					// Map Gambit element types to our types
					var etype mesh.ElementType
					switch gambitType {
					case 1: // Line (1D) - 2 nodes
						etype = mesh.Line
					case 2: // Bar (1D)
						etype = mesh.Line
					case 3: // Triangle (2D)
						etype = mesh.Triangle
					case 4: // Brick (Hex)
						etype = mesh.Hex
					case 5: // Wedge (Prism)
						etype = mesh.Prism
					case 6: // Tetrahedron
						etype = mesh.Tet
					case 7: // Pyramid
						etype = mesh.Pyramid
					default:
						// Skip unknown element types
						continue
					}

					// Read node connectivity
					if len(fields) >= 3+numNodes {
						nodes := make([]int, numNodes)
						for j := 0; j < numNodes; j++ {
							nodeID, _ := strconv.Atoi(fields[3+j])
							// Convert from 1-based to 0-based
							nodes[j] = nodeID - 1
						}

						msh.EtoV = append(msh.EtoV, nodes)
						msh.ElementTypes = append(msh.ElementTypes, etype)
						msh.ElementTags = append(msh.ElementTags, []int{0}) // Default tag
					}
				}
			}

		} else if strings.Contains(line, "ELEMENT GROUP") {
			// Read element groups
			for groupIdx := 0; groupIdx < ngrps; groupIdx++ {
				// Read GROUP line
				if !scanner.Scan() {
					break
				}
				groupLine := strings.TrimSpace(scanner.Text())
				if !strings.HasPrefix(groupLine, "GROUP:") {
					// Try to find next GROUP line
					for scanner.Scan() {
						groupLine = strings.TrimSpace(scanner.Text())
						if strings.HasPrefix(groupLine, "GROUP:") {
							break
						}
						if groupLine == "ENDOFSECTION" {
							break
						}
					}
				}

				if groupLine == "ENDOFSECTION" {
					break
				}

				var groupID, numElems, materialID, nflags int
				parts := strings.Fields(groupLine)

				for i := 0; i < len(parts)-1; i++ {
					if parts[i] == "GROUP:" && i+1 < len(parts) {
						groupID, _ = strconv.Atoi(parts[i+1])
					}
					if parts[i] == "ELEMENTS:" && i+1 < len(parts) {
						numElems, _ = strconv.Atoi(parts[i+1])
					}
					if parts[i] == "MATERIAL:" && i+1 < len(parts) {
						materialID, _ = strconv.Atoi(parts[i+1])
					}
					if parts[i] == "NFLAGS:" && i+1 < len(parts) {
						nflags, _ = strconv.Atoi(parts[i+1])
					}
				}

				// Read entity name
				if !scanner.Scan() {
					break
				}
				entityName := strings.TrimSpace(scanner.Text())

				// Read flags if present
				var flags []int
				if nflags > 0 {
					if !scanner.Scan() {
						break
					}
					flagLine := strings.Fields(scanner.Text())
					flags = make([]int, len(flagLine))
					for i, f := range flagLine {
						flags[i], _ = strconv.Atoi(f)
					}
				}

				// Create element group
				group := &mesh.ElementGroup{
					Dimension:  3, // Gambit groups are typically 3D
					Tag:        groupID,
					Name:       entityName,
					MaterialID: materialID,
					Flags:      flags,
					Elements:   []int{},
				}
				msh.ElementGroups[groupID] = group

				// Read element IDs in this group
				elementsRead := 0
				for elementsRead < numElems {
					if !scanner.Scan() {
						break
					}
					line = strings.TrimSpace(scanner.Text())
					if line == "ENDOFSECTION" {
						break
					}

					fields := strings.Fields(line)
					for _, field := range fields {
						elemID, err := strconv.Atoi(field)
						if err == nil && elemID > 0 && elemID <= len(msh.EtoV) {
							// Elements are 1-indexed in file, 0-indexed in mesh
							elemIdx := elemID - 1
							msh.ElementTags[elemIdx] = []int{groupID}
							group.Elements = append(group.Elements, elemIdx)
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
			for bcIdx := 0; bcIdx < nbsets; bcIdx++ {
				if !scanner.Scan() {
					break
				}
				bcLine := strings.TrimSpace(scanner.Text())

				// Parse boundary condition line
				// Format: NAME ITYPE NENTRY NVALUES IBCODE1 IBCODE2 IBCODE3 IBCODE4 IBCODE5
				parts := strings.Fields(bcLine)
				if len(parts) < 4 {
					continue
				}

				bcName := parts[0]
				itype, _ := strconv.Atoi(parts[1])  // 0=node, 1=element/cell
				nentry, _ := strconv.Atoi(parts[2]) // Number of entries
				// nvalues, _ := strconv.Atoi(parts[3]) // Values per entry

				// Store boundary condition name
				msh.BoundaryTags[bcIdx] = bcName

				// Read boundary data
				if itype == 1 {
					// Element/face boundary conditions
					for i := 0; i < nentry; i++ {
						if !scanner.Scan() {
							break
						}

						fields := strings.Fields(scanner.Text())
						if len(fields) >= 3 {
							elemID, _ := strconv.Atoi(fields[0])
							elemType, _ := strconv.Atoi(fields[1])
							_ = elemType
							faceID, _ := strconv.Atoi(fields[2])

							// Convert element ID from 1-based to 0-based
							elemIdx := elemID - 1

							// Determine boundary element type based on parent element
							if elemIdx >= 0 && elemIdx < len(msh.ElementTypes) {
								parentType := msh.ElementTypes[elemIdx]
								var boundaryType mesh.ElementType

								// Map face to boundary element type
								switch parentType {
								case mesh.Tet:
									boundaryType = mesh.Triangle
								case mesh.Hex:
									boundaryType = mesh.Quad
								case mesh.Prism:
									if faceID <= 2 {
										boundaryType = mesh.Triangle
									} else {
										boundaryType = mesh.Quad
									}
								case mesh.Pyramid:
									if faceID == 1 {
										boundaryType = mesh.Quad
									} else {
										boundaryType = mesh.Triangle
									}
								}

								// Get face vertices
								faceVerts := mesh.GetElementFaces(parentType, msh.EtoV[elemIdx])
								if faceID > 0 && faceID <= len(faceVerts) {
									belem := mesh.BoundaryElement{
										ElementType:   boundaryType,
										Nodes:         faceVerts[faceID-1], // Face IDs are 1-based
										ParentElement: elemIdx,
										ParentFace:    faceID - 1, // Store as 0-based
									}
									msh.AddBoundaryElement(bcName, belem)
								}
							}
						}
					}
				} else {
					// Node boundary conditions - skip reading for now
					for i := 0; i < nentry; i++ {
						if !scanner.Scan() {
							break
						}
					}
				}
			}
		}
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading file: %v", err)
	}

	msh.NumElements = len(msh.EtoV)
	msh.NumVertices = len(msh.Vertices)
	msh.BuildConnectivity()

	return msh, nil
}

package readers

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"

	"github.com/notargets/gocfd/DG3D/mesh"
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

					// Store in 0-based array
					idx := i
					msh.Vertices[idx] = []float64{x, y, z}
					msh.NodeIDMap[nodeID] = idx
					msh.NodeArrayMap[idx] = nodeID
				}
			}

		} else if strings.Contains(line, "ELEMENTS/CELLS") {
			// Read elements
			msh.EtoV = make([][]int, nelem)
			msh.ElementTypes = make([]mesh.ElementType, nelem)
			msh.ElementTags = make([][]int, nelem)
			msh.ElementIDMap = make(map[int]int)

			for i := 0; i < nelem; i++ {
				if !scanner.Scan() {
					return nil, fmt.Errorf("unexpected EOF reading elements")
				}

				fields := strings.Fields(scanner.Text())
				if len(fields) >= 3 {
					elemID, _ := strconv.Atoi(fields[0])
					elemType, _ := strconv.Atoi(fields[1])
					numNodes, _ := strconv.Atoi(fields[2])

					// Map Gambit element types to mesh types
					var meshType mesh.ElementType
					switch elemType {
					case 4: // Hexahedron
						meshType = mesh.Hex
					case 5: // Prism/Wedge
						meshType = mesh.Prism
					case 6: // Tetrahedron
						meshType = mesh.Tet
					case 7: // Pyramid
						meshType = mesh.Pyramid
					case 3: // Quadrilateral (2D)
						meshType = mesh.Quad
					case 2: // Triangle (2D)
						meshType = mesh.Triangle
					default:
						continue // Skip unsupported element types
					}

					msh.ElementTypes[i] = meshType
					msh.ElementIDMap[elemID] = i

					// Read node connectivity
					nodes := make([]int, 0, numNodes)
					for j := 3; j < len(fields) && len(nodes) < numNodes; j++ {
						nodeID, err := strconv.Atoi(fields[j])
						if err == nil {
							// Convert node ID to 0-based index
							if idx, ok := msh.NodeIDMap[nodeID]; ok {
								nodes = append(nodes, idx)
							}
						}
					}
					msh.EtoV[i] = nodes

					// Initialize with default tag
					msh.ElementTags[i] = []int{0}
				}
			}

		} else if strings.Contains(line, "ELEMENT GROUP") {
			// Initialize ElementGroups map if needed
			if msh.ElementGroups == nil {
				msh.ElementGroups = make(map[int]*mesh.ElementGroup)
			}

			// Read element group information
			for grpIdx := 0; grpIdx < ngrps; grpIdx++ {
				if !scanner.Scan() {
					break
				}
				groupLine := strings.TrimSpace(scanner.Text())

				// Parse GROUP line
				// Format: GROUP: NGP ELEMENTS: NELGP MATERIAL: MTYP NFLAGS: NFLAGS
				var groupID, numElems, material, flags int
				parts := strings.Fields(groupLine)
				for i := 0; i < len(parts)-1; i++ {
					if parts[i] == "GROUP:" && i+1 < len(parts) {
						groupID, _ = strconv.Atoi(parts[i+1])
					} else if parts[i] == "ELEMENTS:" && i+1 < len(parts) {
						numElems, _ = strconv.Atoi(parts[i+1])
					} else if parts[i] == "MATERIAL:" && i+1 < len(parts) {
						material, _ = strconv.Atoi(parts[i+1])
					} else if parts[i] == "NFLAGS:" && i+1 < len(parts) {
						flags, _ = strconv.Atoi(parts[i+1])
					}
				}

				_ = material
				// Read entity name
				var entityName string
				if scanner.Scan() {
					entityName = strings.TrimSpace(scanner.Text())
				}

				// Skip flags if present
				if flags > 0 && scanner.Scan() {
					// Skip flags line
				}

				// Create element group
				group := &mesh.ElementGroup{
					Tag:      groupID,
					Name:     entityName,
					Elements: []int{},
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
						if err == nil {
							// Convert element ID to 0-based index
							if elemIdx, ok := msh.ElementIDMap[elemID]; ok {
								msh.ElementTags[elemIdx] = []int{groupID}
								group.Elements = append(group.Elements, elemIdx)
							}
						}
						elementsRead++
						if elementsRead >= numElems {
							break
						}
					}
				}
			}

		} else if strings.Contains(line, "BOUNDARY CONDITIONS") {
			// Initialize BoundaryTags map if needed
			if msh.BoundaryTags == nil {
				msh.BoundaryTags = make(map[int]string)
			}

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
					return nil, fmt.Errorf("invalid boundary condition format at line %d: expected at least 4 fields, got %d", bcIdx+1, len(parts))
				}

				// Validate that fields 2-4 are integers
				bcName := parts[0]
				itype, err1 := strconv.Atoi(parts[1])  // 0=node, 1=element/cell
				nentry, err2 := strconv.Atoi(parts[2]) // Number of entries
				_, err3 := strconv.Atoi(parts[3])      // Values per entry (nvalues)

				if err1 != nil || err2 != nil || err3 != nil {
					// Provide more specific error messages
					if err2 != nil {
						return nil, fmt.Errorf("invalid boundary condition format: NENTRY must be an integer, got '%s'", parts[2])
					}
					return nil, fmt.Errorf("invalid boundary condition format: ITYPE, NENTRY, and NVALUES must be integers")
				}

				// Validate ITYPE is 0 or 1
				if itype != 0 && itype != 1 {
					return nil, fmt.Errorf("invalid boundary condition ITYPE: %d (must be 0 for node or 1 for element/cell)", itype)
				}

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
							elemID, err1 := strconv.Atoi(fields[0])
							elemType, err2 := strconv.Atoi(fields[1])
							faceID, err3 := strconv.Atoi(fields[2])

							if err1 != nil || err2 != nil || err3 != nil {
								return nil, fmt.Errorf("invalid element boundary data: expected integers for element ID, type, and face ID")
							}

							_ = elemType

							// Convert element ID to 0-based index
							if elemIdx, ok := msh.ElementIDMap[elemID]; ok {
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

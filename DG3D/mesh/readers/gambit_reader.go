package readers

import (
	"bufio"
	"fmt"
	"github.com/notargets/gocfd/utils"
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
			msh.ElementTypes = make([]utils.ElementType, nelem)
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
					var meshType utils.ElementType
					switch elemType {
					case 4: // Hexahedron
						meshType = utils.Hex
					case 5: // Prism/Wedge
						meshType = utils.Prism
					case 6: // Tetrahedron
						meshType = utils.Tet
					case 7: // Pyramid
						meshType = utils.Pyramid
					case 3: // Quadrilateral (2D)
						meshType = utils.Quad
					case 2: // Triangle (2D)
						meshType = utils.Triangle
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

			// Initialize EToP array if we have partitions (ngrps > 1)
			if ngrps > 1 && msh.EToP == nil {
				msh.EToP = make([]int, nelem)
				// Initialize all elements to partition -1 (unassigned)
				for i := range msh.EToP {
					msh.EToP[i] = -1
				}
			}

			// Read GROUP line
			if !scanner.Scan() {
				continue
			}
			groupLine := scanner.Text()

			// Parse GROUP line: GROUP: <id> ELEMENTS: <count> MATERIAL: <mat> NFLAGS: <flags>
			var groupID, elemCount int
			fields := strings.Fields(groupLine)
			for i := 0; i < len(fields); i++ {
				if fields[i] == "GROUP:" && i+1 < len(fields) {
					groupID, _ = strconv.Atoi(fields[i+1])
				} else if fields[i] == "ELEMENTS:" && i+1 < len(fields) {
					elemCount, _ = strconv.Atoi(fields[i+1])
				}
			}

			// Read group name
			if !scanner.Scan() {
				continue
			}
			groupName := strings.TrimSpace(scanner.Text())

			// Skip flags line
			if !scanner.Scan() {
				continue
			}

			// Create element group
			group := &mesh.ElementGroup{
				Tag:       groupID,
				Name:      groupName,
				Dimension: 3, // Assuming 3D for now
				Elements:  make([]int, 0, elemCount),
			}

			// Determine partition ID from group name if it matches pattern p0, p1, p2, etc.
			partitionID := -1
			if strings.HasPrefix(groupName, "p") && len(groupName) > 1 {
				// Try to parse partition number
				partNum, err := strconv.Atoi(groupName[1:])
				if err == nil {
					partitionID = partNum
				}
			}

			// Read element IDs
			elementsRead := 0
			for scanner.Scan() && elementsRead < elemCount {
				line := strings.TrimSpace(scanner.Text())
				if line == "ENDOFSECTION" {
					break
				}

				fields := strings.Fields(line)
				for _, field := range fields {
					elemID, err := strconv.Atoi(field)
					if err == nil && elemID > 0 {
						// Convert to 0-based index
						if idx, ok := msh.ElementIDMap[elemID]; ok {
							group.Elements = append(group.Elements, idx)

							// If this is a partition group, update EToP
							if partitionID >= 0 && msh.EToP != nil {
								msh.EToP[idx] = partitionID
							}
						}
						elementsRead++
						if elementsRead >= elemCount {
							break
						}
					}
				}
			}

			msh.ElementGroups[groupID] = group

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
								return nil, fmt.Errorf("invalid boundary condition format")
							}

							// Convert element ID to 0-based index
							if idx, ok := msh.ElementIDMap[elemID]; ok {
								// Initialize BoundaryElements if needed
								if msh.BoundaryElements == nil {
									msh.BoundaryElements = make(map[string][]mesh.BoundaryElement)
								}

								// Create boundary element
								be := mesh.BoundaryElement{
									ParentElement: idx,
									ParentFace:    faceID - 1, // Convert to 0-based
								}

								// Map element type to boundary element type
								switch elemType {
								case 3:
									be.ElementType = utils.Quad
								case 2:
									be.ElementType = utils.Triangle
								case 10:
									be.ElementType = utils.Triangle
								case 11:
									be.ElementType = utils.Quad
								}

								msh.BoundaryElements[bcName] = append(msh.BoundaryElements[bcName], be)
							}
						}
					}
				} else {
					// Node boundary conditions - skip for now
					for i := 0; i < nentry; i++ {
						if !scanner.Scan() {
							break
						}
					}
				}
			}
		}
	}

	// Set mesh properties
	msh.NumVertices = len(msh.Vertices)
	msh.NumElements = len(msh.EtoV)

	return msh, scanner.Err()
}

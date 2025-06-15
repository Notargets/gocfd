package readers

import (
	"bufio"
	"fmt"
	"github.com/notargets/gocfd/DG3D/mesh"
	"os"
	"strconv"
	"strings"
)

// ReadGambitNeutral reads a Gambit neutral file
func ReadGambitNeutral(filename string) (*mesh.Mesh, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	msh := mesh.NewMesh()
	scanner := bufio.NewScanner(file)

	var numnp, nelem, ngrps, nbsets int

	// Read until we find the problem size parameters
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())

		// Look for the problem size header line
		if strings.Contains(line, "NUMNP") && strings.Contains(line, "NELEM") {
			// Next line contains the actual values
			if scanner.Scan() {
				values := strings.Fields(scanner.Text())
				if len(values) >= 4 {
					numnp, _ = strconv.Atoi(values[0])
					nelem, _ = strconv.Atoi(values[1])
					ngrps, _ = strconv.Atoi(values[2])
					nbsets, _ = strconv.Atoi(values[3])
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
			msh.Vertices = make([][]float64, numnp)
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
						msh.Vertices[id-1] = []float64{x, y, z}
						nodeCount++
					}
				}
			}

		} else if strings.Contains(line, "ELEMENTS/CELLS") {
			// Read elements
			msh.EtoV = make([][]int, 0, nelem)
			msh.ElementTypes = make([]mesh.ElementType, 0, nelem)
			msh.ElementTags = make([][]int, 0, nelem)

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
					var etype mesh.ElementType
					validElem := true
					switch elemType {
					case 1: // Edge - skip
						validElem = false
					case 2: // Quad - skip 2D
						validElem = false
					case 3: // Triangle - skip 2D
						validElem = false
					case 4: // Brick
						etype = mesh.Hex
					case 5: // Wedge/Prism
						etype = mesh.Prism
					case 6: // Tetrahedron
						etype = mesh.Tet
					case 7: // Pyramid
						etype = mesh.Pyramid
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

						msh.EtoV = append(msh.EtoV, verts)
						msh.ElementTypes = append(msh.ElementTypes, etype)
						msh.ElementTags = append(msh.ElementTags, []int{0}) // Default tag
					}
				}
			}

		} else if strings.Contains(line, "ELEMENT GROUP") {
			// Read element groups - process all groups
			for groupIdx := 0; groupIdx < ngrps; groupIdx++ {
				// Read the group header line
				if !scanner.Scan() {
					break
				}
				groupLine := strings.TrimSpace(scanner.Text())

				// Parse GROUP: NGP ELEMENTS: NELGP MATERIAL: MTYP NFLAGS: NFLAGS
				var groupID, numElems, nflags int
				parts := strings.Fields(groupLine)

				for i := 0; i < len(parts)-1; i++ {
					if parts[i] == "GROUP:" && i+1 < len(parts) {
						groupID, _ = strconv.Atoi(parts[i+1])
					}
					if parts[i] == "ELEMENTS:" && i+1 < len(parts) {
						numElems, _ = strconv.Atoi(parts[i+1])
					}
					if parts[i] == "NFLAGS:" && i+1 < len(parts) {
						nflags, _ = strconv.Atoi(parts[i+1])
					}
				}

				// Read entity name
				if !scanner.Scan() {
					break
				}
				// entityName := strings.TrimSpace(scanner.Text())

				// Read flags if present
				if nflags > 0 {
					if !scanner.Scan() {
						break
					}
					// flags := strings.TrimSpace(scanner.Text())
				}

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
							msh.ElementTags[elemID-1] = []int{groupID}
						}
						elementsRead++
						if elementsRead >= numElems {
							break
						}
					}
				}

				// Check if we've reached ENDOFSECTION
				if line == "ENDOFSECTION" {
					break
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
				parts := strings.Fields(bcLine)
				if len(parts) >= 2 {
					bcName := parts[0]
					// Could parse ITYPE, NENTRY, etc. if needed

					// Store boundary condition name
					msh.BoundaryTags[bcIdx] = bcName
				}

				// Skip to ENDOFSECTION or next BC
				// In a full implementation, we would read the BC data here
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

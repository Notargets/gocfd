package readers

import (
	"bufio"
	"fmt"
	"github.com/notargets/gocfd/DG3D/mesh"
	"github.com/notargets/gocfd/utils"
	"os"
	"strconv"
	"strings"
)

// ReadSU2 reads an SU2 native format file
func ReadSU2(filename string) (*mesh.Mesh, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	msh := mesh.NewMesh()
	scanner := bufio.NewScanner(file)

	var ndime int
	var hasNDIME, hasNPOIN bool

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())

		// Skip empty lines
		if line == "" {
			continue
		}

		// Skip comments (text after %)
		if idx := strings.Index(line, "%"); idx >= 0 {
			line = strings.TrimSpace(line[:idx])
			if line == "" {
				continue
			}
		}

		if strings.HasPrefix(line, "NDIME=") {
			hasNDIME = true
			fmt.Sscanf(line, "NDIME=%d", &ndime)
			if ndime != 2 && ndime != 3 {
				return nil, fmt.Errorf("unsupported dimension: NDIME=%d", ndime)
			}

		} else if strings.HasPrefix(line, "NPOIN=") {
			hasNPOIN = true
			var npoin int
			fmt.Sscanf(line, "NPOIN=%d", &npoin)

			msh.Vertices = make([][]float64, npoin)
			msh.NodeIDMap = make(map[int]int)
			msh.NodeArrayMap = make(map[int]int)

			for i := 0; i < npoin; i++ {
				if !scanner.Scan() {
					return nil, fmt.Errorf("unexpected EOF reading nodes")
				}

				fields := strings.Fields(scanner.Text())
				if len(fields) < ndime {
					return nil, fmt.Errorf("invalid node line: expected at least %d coordinates", ndime)
				}

				coords := make([]float64, 3) // Always store 3D coordinates
				for j := 0; j < ndime; j++ {
					coords[j], err = strconv.ParseFloat(fields[j], 64)
					if err != nil {
						return nil, fmt.Errorf("invalid coordinate: %v", err)
					}
				}

				// Node ID is implicit (0-based) based on order
				// Legacy format may have explicit ID at end of line (ignored)
				msh.Vertices[i] = coords
				msh.NodeIDMap[i] = i
				msh.NodeArrayMap[i] = i
			}

		} else if strings.HasPrefix(line, "NELEM=") {
			var nelem int
			fmt.Sscanf(line, "NELEM=%d", &nelem)

			// Pre-allocate slices
			msh.EtoV = make([][]int, 0, nelem)
			msh.ElementTypes = make([]utils.ElementType, 0, nelem)
			msh.ElementTags = make([][]int, 0, nelem)

			// Read elements
			for i := 0; i < nelem; i++ {
				if !scanner.Scan() {
					return nil, fmt.Errorf("unexpected EOF reading elements")
				}

				fields := strings.Fields(scanner.Text())
				if len(fields) < 2 {
					return nil, fmt.Errorf("invalid element line")
				}

				su2Type, err := strconv.Atoi(fields[0])
				if err != nil {
					return nil, fmt.Errorf("invalid element type: %v", err)
				}

				// Map SU2/VTK element types to our types
				etype, ok := su2ElementTypeMap[su2Type]
				if !ok {
					return nil, fmt.Errorf("unknown element type: %d", su2Type)
				}

				numNodes := etype.GetNumNodes()
				if len(fields) < numNodes+1 {
					return nil, fmt.Errorf("element type %v expects %d nodes, got %d fields",
						etype, numNodes, len(fields)-1)
				}

				// Read node indices
				nodes := make([]int, numNodes)
				for j := 0; j < numNodes; j++ {
					nodes[j], err = strconv.Atoi(fields[1+j])
					if err != nil {
						return nil, fmt.Errorf("invalid node index: %v", err)
					}
					// Validate node index
					if nodes[j] < 0 || nodes[j] >= len(msh.Vertices) {
						return nil, fmt.Errorf("node index %d out of range [0,%d)",
							nodes[j], len(msh.Vertices))
					}
				}

				// Element ID is implicit (0-based) based on order
				// Legacy format may have explicit ID at end of line (ignored)
				msh.EtoV = append(msh.EtoV, nodes)
				msh.ElementTypes = append(msh.ElementTypes, etype)
				msh.ElementTags = append(msh.ElementTags, []int{0}) // Default tag
			}

		} else if strings.HasPrefix(line, "NMARK=") {
			var nmark int
			fmt.Sscanf(line, "NMARK=%d", &nmark)

			msh.BoundaryTags = make(map[int]string)

			// Read boundary markers
			for i := 0; i < nmark; i++ {
				// Read MARKER_TAG= line
				if !scanner.Scan() {
					return nil, fmt.Errorf("unexpected EOF reading marker %d", i)
				}

				markerLine := strings.TrimSpace(scanner.Text())
				if !strings.HasPrefix(markerLine, "MARKER_TAG=") {
					return nil, fmt.Errorf("expected MARKER_TAG=, got: %s", markerLine)
				}

				tagName := strings.TrimSpace(strings.TrimPrefix(markerLine, "MARKER_TAG="))

				// Read MARKER_ELEMS= line
				if !scanner.Scan() {
					return nil, fmt.Errorf("unexpected EOF reading marker elements for %s", tagName)
				}

				elemLine := strings.TrimSpace(scanner.Text())
				var nMarkerElems int
				if _, err := fmt.Sscanf(elemLine, "MARKER_ELEMS=%d", &nMarkerElems); err != nil {
					return nil, fmt.Errorf("invalid MARKER_ELEMS line: %s", elemLine)
				}

				// Store boundary tag
				msh.BoundaryTags[i] = tagName

				// Read boundary elements
				for j := 0; j < nMarkerElems; j++ {
					if !scanner.Scan() {
						return nil, fmt.Errorf("unexpected EOF reading boundary elements")
					}

					fields := strings.Fields(scanner.Text())
					if len(fields) < 2 {
						return nil, fmt.Errorf("invalid boundary element line")
					}

					// Parse boundary element type
					boundaryType, err := strconv.Atoi(fields[0])
					if err != nil {
						return nil, fmt.Errorf("invalid boundary element type: %v", err)
					}

					// Map boundary element type
					var btype utils.ElementType
					switch boundaryType {
					case 3: // Line (2D boundary)
						btype = utils.Line
					case 5: // Triangle (3D boundary)
						btype = utils.Triangle
					case 9: // Quad (3D boundary)
						btype = utils.Quad
					default:
						return nil, fmt.Errorf("unknown boundary element type: %d", boundaryType)
					}

					numBoundaryNodes := btype.GetNumNodes()
					if len(fields) < numBoundaryNodes+1 {
						return nil, fmt.Errorf("boundary element type %v expects %d nodes",
							btype, numBoundaryNodes)
					}

					// Read boundary nodes
					boundaryNodes := make([]int, numBoundaryNodes)
					for k := 0; k < numBoundaryNodes; k++ {
						boundaryNodes[k], err = strconv.Atoi(fields[1+k])
						if err != nil {
							return nil, fmt.Errorf("invalid boundary node index: %v", err)
						}
					}

					// Create boundary element
					boundaryElem := mesh.BoundaryElement{
						ElementType:   btype,
						Nodes:         boundaryNodes,
						ParentElement: -1, // Not tracked in SU2 format
						ParentFace:    -1, // Not tracked in SU2 format
					}

					msh.AddBoundaryElement(tagName, boundaryElem)
				}
			}
		}
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading file: %v", err)
	}

	// Validate that we read the required sections
	if !hasNDIME {
		return nil, fmt.Errorf("missing required NDIME= section")
	}

	if !hasNPOIN {
		return nil, fmt.Errorf("missing required NPOIN= section")
	}

	msh.NumElements = len(msh.EtoV)
	msh.NumVertices = len(msh.Vertices)
	msh.BuildConnectivity()

	return msh, nil
}

// su2ElementTypeMap maps SU2/VTK element type identifiers to our ElementType
var su2ElementTypeMap = map[int]utils.ElementType{
	3:  utils.Line,     // VTK_LINE
	5:  utils.Triangle, // VTK_TRIANGLE
	9:  utils.Quad,     // VTK_QUAD
	10: utils.Tet,      // VTK_TETRA
	12: utils.Hex,      // VTK_HEXAHEDRON
	13: utils.Prism,    // VTK_WEDGE
	14: utils.Pyramid,  // VTK_PYRAMID
}

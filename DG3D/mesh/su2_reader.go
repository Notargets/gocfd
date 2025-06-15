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

			mesh.Vertices = make([][]float64, npoin)
			mesh.NodeIDMap = make(map[int]int)
			mesh.NodeArrayMap = make(map[int]int)

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
				mesh.Vertices[i] = coords
				mesh.NodeIDMap[i] = i
				mesh.NodeArrayMap[i] = i
			}

		} else if strings.HasPrefix(line, "NELEM=") {
			var nelem int
			fmt.Sscanf(line, "NELEM=%d", &nelem)

			// Pre-allocate slices
			mesh.EtoV = make([][]int, 0, nelem)
			mesh.ElementTypes = make([]ElementType, 0, nelem)
			mesh.ElementTags = make([][]int, 0, nelem)

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
					if nodes[j] < 0 || nodes[j] >= len(mesh.Vertices) {
						return nil, fmt.Errorf("node index %d out of range [0,%d)",
							nodes[j], len(mesh.Vertices))
					}
				}

				// Element ID is implicit (0-based) based on order
				// Legacy format may have explicit ID at end of line (ignored)
				mesh.EtoV = append(mesh.EtoV, nodes)
				mesh.ElementTypes = append(mesh.ElementTypes, etype)
				mesh.ElementTags = append(mesh.ElementTags, []int{0}) // Default tag
			}

		} else if strings.HasPrefix(line, "NMARK=") {
			var nmark int
			fmt.Sscanf(line, "NMARK=%d", &nmark)

			mesh.BoundaryTags = make(map[int]string)

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
				mesh.BoundaryTags[i] = tagName

				// Read boundary elements
				// For now, we skip the actual boundary element data
				// A full implementation would store these in a separate structure
				for j := 0; j < nMarkerElems; j++ {
					if !scanner.Scan() {
						return nil, fmt.Errorf("unexpected EOF reading boundary elements")
					}
					// Parse boundary element if needed
					// fields := strings.Fields(scanner.Text())
					// boundaryType := fields[0]
					// boundaryNodes := fields[1:]
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

	mesh.NumElements = len(mesh.EtoV)
	mesh.NumVertices = len(mesh.Vertices)
	mesh.BuildConnectivity()

	return mesh, nil
}

// su2ElementTypeMap maps SU2/VTK element type identifiers to our ElementType
var su2ElementTypeMap = map[int]ElementType{
	3:  Line,     // VTK_LINE
	5:  Triangle, // VTK_TRIANGLE
	9:  Quad,     // VTK_QUAD
	10: Tet,      // VTK_TETRA
	12: Hex,      // VTK_HEXAHEDRON
	13: Prism,    // VTK_WEDGE
	14: Pyramid,  // VTK_PYRAMID
}

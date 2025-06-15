package readers

import (
	"bufio"
	"fmt"
	"github.com/notargets/gocfd/DG3D/mesh"
	"os"
	"strings"
)

// ReadGmshAuto automatically detects the Gmsh format version and reads the file
func ReadGmshAuto(filename string) (*mesh.Mesh, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	var version string

	// Look for $MeshFormat section to determine version
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "$MeshFormat" {
			if scanner.Scan() {
				parts := strings.Fields(scanner.Text())
				if len(parts) > 0 {
					version = parts[0]
					break
				}
			}
		}
	}

	file.Close()

	// Determine which reader to use based on version
	if strings.HasPrefix(version, "4.") {
		return ReadGmsh4(filename)
	} else if strings.HasPrefix(version, "2.") {
		return ReadGmsh22(filename)
	} else if version == "" {
		return nil, fmt.Errorf("could not find $MeshFormat section")
	} else {
		return nil, fmt.Errorf("unsupported Gmsh format version: %s", version)
	}
}

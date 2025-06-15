package readers

import (
	"fmt"
	"github.com/notargets/gocfd/DG3D/mesh"
	"path/filepath"
	"strings"
)

// ReadMeshFile reads a mesh file based on extension
func ReadMeshFile(filename string) (*mesh.Mesh, error) {
	ext := strings.ToLower(filepath.Ext(filename))

	switch ext {
	case ".neu":
		return ReadGambitNeutral(filename)
	case ".msh":
		// Try to detect version by reading first few lines
		return ReadGmshAuto(filename)
	case ".su2":
		return ReadSU2(filename)
	default:
		return nil, fmt.Errorf("unsupported mesh format: %s", ext)
	}
}

package readfiles

import (
	"fmt"

	"github.com/notargets/gocfd/geometry2D"
)

func PrintTriMesh(trimesh *geometry2D.TriMesh) {
	fmt.Println("Geometry, X, Y for each vertex:")
	for _, pt := range trimesh.Geometry {
		fmt.Printf("%5.2f %5.2f\n", pt.X[0], pt.X[1])
	}
	fmt.Println("Vertex indices, 3 for each triangle:")
	for _, t := range trimesh.Triangles {
		for i := 0; i < 3; i++ {
			fmt.Printf("%d ", t.Nodes[i])
		}
		fmt.Printf("\n")
	}
	fmt.Println("Scalar color value, one for each index in the triangle:")
	for _, a := range trimesh.Attributes {
		for i := 0; i < 3; i++ {
			fmt.Printf("%5.2f ", a[i])
		}
		fmt.Printf("\n")
	}
}

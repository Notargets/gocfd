package main

import (
	"flag"
	"fmt"
	"log"
	"os"

	"github.com/notargets/gocfd/DG3D/mesh"
)

func main() {
	// Command line flags
	var (
		meshFile  = flag.String("mesh", "", "Input mesh file (.neu, .msh, or .su2)")
		nparts    = flag.Int("nparts", 4, "Number of partitions")
		imbalance = flag.Float64("imbalance", 0.05, "Allowed load imbalance (0.05 = 5%)")
		objective = flag.String("obj", "vol", "Objective: 'cut' or 'vol'")
		output    = flag.String("output", "", "Output file for partitioned mesh")
		verbose   = flag.Bool("v", false, "Verbose output")
	)
	flag.Parse()

	if *meshFile == "" {
		fmt.Fprintf(os.Stderr, "Error: mesh file required\n")
		flag.Usage()
		os.Exit(1)
	}

	// Read mesh - use different variable name to avoid shadowing package
	log.Printf("Reading mesh from %s", *meshFile)
	m, err := mesh.ReadMeshFile(*meshFile)
	if err != nil {
		log.Fatalf("Failed to read mesh: %v", err)
	}

	// Print mesh statistics
	if *verbose {
		m.PrintStatistics()
	}

	// Create partition configuration
	config := &mesh.PartitionConfig{
		NumPartitions:    int32(*nparts),
		ImbalanceFactor:  float32(1.0 + *imbalance),
		UseEdgeWeights:   true,
		UseVertexWeights: true,
		Objective:        *objective,
	}

	// Create partitioner
	partitioner := mesh.NewMeshPartitioner(m, config)

	// Perform partitioning
	log.Printf("Partitioning mesh into %d parts", *nparts)
	err = partitioner.Partition()
	if err != nil {
		log.Fatalf("Partitioning failed: %v", err)
	}

	// Print partition assignment if verbose
	if *verbose {
		fmt.Println("\nElement partition assignment:")
		for i := 0; i < m.NumElements; i++ {
			fmt.Printf("  Element %d (type=%s) -> Partition %d\n",
				i, m.ElementTypes[i], m.EToP[i])
		}
	}

	// Get partition boundary information
	boundaryFaces := partitioner.GetPartitionBoundaryFaces()
	fmt.Printf("\nPartition boundary faces:\n")
	for part, faces := range boundaryFaces {
		fmt.Printf("  Partition %d: %d boundary faces\n", part, len(faces))
	}

	// Export if requested
	if *output != "" {
		log.Printf("Exporting partitioned mesh to %s", *output)
		err = exportPartitionedMesh(m, *output)
		if err != nil {
			log.Printf("Warning: Failed to export mesh: %v", err)
		}
	}

	log.Println("Partitioning complete!")
}

// exportPartitionedMesh writes a simple format with partition information
func exportPartitionedMesh(m *mesh.Mesh, filename string) error {
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	// Write header
	fmt.Fprintf(file, "# Partitioned Mesh\n")
	fmt.Fprintf(file, "# Vertices: %d\n", m.NumVertices)
	fmt.Fprintf(file, "# Elements: %d\n", m.NumElements)
	fmt.Fprintf(file, "# Partitions: %d\n", getMaxPartition(m.EToP)+1)
	fmt.Fprintf(file, "\n")

	// Write vertices
	fmt.Fprintf(file, "VERTICES %d\n", m.NumVertices)
	for i, v := range m.Vertices {
		fmt.Fprintf(file, "%d %.6f %.6f %.6f\n", i, v[0], v[1], v[2])
	}
	fmt.Fprintf(file, "\n")

	// Write elements with partition assignment
	fmt.Fprintf(file, "ELEMENTS %d\n", m.NumElements)
	for i := 0; i < m.NumElements; i++ {
		fmt.Fprintf(file, "%d %s %d", i, m.ElementTypes[i], m.EToP[i])
		for _, v := range m.Elements[i] {
			fmt.Fprintf(file, " %d", v)
		}
		fmt.Fprintf(file, "\n")
	}

	return nil
}

func getMaxPartition(EToP []int) int {
	max := 0
	for _, p := range EToP {
		if p > max {
			max = p
		}
	}
	return max
}

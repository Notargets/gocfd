package mesh

import (
	"fmt"
	"log"
	"math"

	metis "github.com/notargets/go-metis"
)

// PartitionConfig holds configuration for mesh partitioning
type PartitionConfig struct {
	NumPartitions    int32
	ImbalanceFactor  float32 // e.g., 1.05 for 5% imbalance
	UseEdgeWeights   bool
	UseVertexWeights bool
	Objective        string // "cut" or "vol"
}

// DefaultPartitionConfig returns default partitioning configuration
func DefaultPartitionConfig(nparts int32) *PartitionConfig {
	return &PartitionConfig{
		NumPartitions:    nparts,
		ImbalanceFactor:  1.05,
		UseEdgeWeights:   true,
		UseVertexWeights: true,
		Objective:        "vol", // minimize communication volume
	}
}

// MeshPartitioner handles partitioning of mixed element meshes
type MeshPartitioner struct {
	mesh   *Mesh
	config *PartitionConfig

	// Cost models
	computeCostModel func(elemType ElementType, numVertices int) int32
	commCostModel    func(faceVertices int, isBoundary bool) int32
}

// NewMeshPartitioner creates a new partition_mesh for the given mesh
func NewMeshPartitioner(mesh *Mesh, config *PartitionConfig) *MeshPartitioner {
	mp := &MeshPartitioner{
		mesh:   mesh,
		config: config,
	}

	// Default compute cost model
	mp.computeCostModel = func(elemType ElementType, numVertices int) int32 {
		// Base costs reflect relative computational expense
		baseCost := map[ElementType]int32{
			Tet:     1,
			Hex:     8, // Hex has 8 vertices vs 4 for tet
			Prism:   6,
			Pyramid: 5,
		}

		// Could be enhanced with polynomial order information
		return baseCost[elemType]
	}

	// Default communication cost model
	mp.commCostModel = func(faceVertices int, isBoundary bool) int32 {
		if isBoundary {
			return 0 // No communication across boundaries
		}

		// Cost proportional to number of face DOFs
		// For linear elements, this is number of vertices
		// For higher order, would be (p+1)^2 for quads, (p+1)(p+2)/2 for triangles
		switch faceVertices {
		case 3: // Triangle face
			return 3
		case 4: // Quad face
			return 4
		default:
			return int32(faceVertices)
		}
	}

	return mp
}

// Partition performs the mesh partitioning
func (mp *MeshPartitioner) Partition() error {
	log.Printf("Partitioning mesh with %d elements into %d parts",
		mp.mesh.NumElements, mp.config.NumPartitions)

	// Build METIS graph
	xadj, adjncy, vwgt, adjwgt := mp.buildMetisGraph()

	// Set METIS options
	opts := make([]int32, metis.NoOptions)
	err := metis.SetDefaultOptions(opts)
	if err != nil {
		return fmt.Errorf("failed to set METIS options: %w", err)
	}

	// Set objective function
	if mp.config.Objective == "vol" {
		opts[metis.OptionObjType] = metis.ObjTypeVol
	} else {
		opts[metis.OptionObjType] = metis.ObjTypeCut
	}

	// Set allowed imbalance
	ubvec := []float32{mp.config.ImbalanceFactor}

	// Handle case where weights might be nil
	var vwgtPtr, adjwgtPtr []int32
	if mp.config.UseVertexWeights {
		vwgtPtr = vwgt
	}
	if mp.config.UseEdgeWeights {
		adjwgtPtr = adjwgt
	}

	// Perform partitioning
	part, objval, err := metis.PartGraphKwayWeighted(
		xadj, adjncy, vwgtPtr, adjwgtPtr,
		mp.config.NumPartitions, nil, ubvec, opts,
	)
	if err != nil {
		return fmt.Errorf("METIS partitioning failed: %w", err)
	}

	// Store partition assignment
	mp.mesh.EToP = make([]int, mp.mesh.NumElements)
	for i := 0; i < mp.mesh.NumElements; i++ {
		mp.mesh.EToP[i] = int(part[i])
	}

	// Analyze partition quality
	mp.analyzePartition(objval)

	return nil
}

// buildMetisGraph converts mesh connectivity to METIS format
func (mp *MeshPartitioner) buildMetisGraph() (xadj, adjncy, vwgt, adjwgt []int32) {
	ne := mp.mesh.NumElements

	// Build vertex weights (computational cost per element)
	if mp.config.UseVertexWeights {
		vwgt = make([]int32, ne)
		for i := 0; i < ne; i++ {
			vwgt[i] = mp.computeCostModel(
				mp.mesh.ElementTypes[i],
				len(mp.mesh.Elements[i]),
			)
		}
	}

	// Build adjacency and edge weights
	xadj = make([]int32, ne+1)
	adjncy = []int32{}
	adjwgt = []int32{}

	xadj[0] = 0
	for elem := 0; elem < ne; elem++ {
		for faceIdx, neighbor := range mp.mesh.EToE[elem] {
			if neighbor >= 0 && neighbor != elem {
				// Add neighbor
				adjncy = append(adjncy, int32(neighbor))

				// Add edge weight (communication cost)
				if mp.config.UseEdgeWeights {
					faceID := mp.mesh.EToF[elem][faceIdx]
					face := mp.mesh.Faces[faceID]
					cost := mp.commCostModel(len(face.Vertices), false)
					adjwgt = append(adjwgt, cost)
				}
			}
		}
		xadj[elem+1] = int32(len(adjncy))
	}

	return xadj, adjncy, vwgt, adjwgt
}

// analyzePartition computes and reports partition quality metrics
func (mp *MeshPartitioner) analyzePartition(objval int32) {
	nparts := int(mp.config.NumPartitions)

	// Initialize partition statistics
	partStats := make([]PartitionStats, nparts)
	for i := range partStats {
		partStats[i].ID = i
		// FIX: Initialize the maps
		partStats[i].ElementTypes = make(map[ElementType]int)
		partStats[i].NumNeighbors = make(map[int]int)
	}

	// Gather element statistics
	for elem := 0; elem < mp.mesh.NumElements; elem++ {
		part := mp.mesh.EToP[elem]
		stats := &partStats[part]

		stats.NumElements++
		stats.ElementTypes[mp.mesh.ElementTypes[elem]]++

		// Add computational cost
		cost := mp.computeCostModel(
			mp.mesh.ElementTypes[elem],
			len(mp.mesh.Elements[elem]),
		)
		stats.ComputeLoad += int64(cost)
	}

	// Analyze communication
	cutEdges := 0
	commVolume := int64(0)
	interfaceFaces := make(map[[2]int][]int) // [part1,part2] -> face IDs

	for elem := 0; elem < mp.mesh.NumElements; elem++ {
		elemPart := mp.mesh.EToP[elem]

		for faceIdx, neighbor := range mp.mesh.EToE[elem] {
			if neighbor >= 0 && neighbor > elem { // Count each edge once
				neighborPart := mp.mesh.EToP[neighbor]

				if elemPart != neighborPart {
					cutEdges++

					// Track interface
					p1, p2 := elemPart, neighborPart
					if p1 > p2 {
						p1, p2 = p2, p1
					}
					faceID := mp.mesh.EToF[elem][faceIdx]
					interfaceFaces[[2]int{p1, p2}] = append(
						interfaceFaces[[2]int{p1, p2}], faceID)

					// Add communication cost
					face := mp.mesh.Faces[faceID]
					cost := mp.commCostModel(len(face.Vertices), false)
					commVolume += int64(cost)

					// Update partition stats
					partStats[elemPart].NumNeighbors[neighborPart]++
					partStats[neighborPart].NumNeighbors[elemPart]++
				}
			}
		}
	}

	// Compute load imbalance
	avgLoad := float64(0)
	maxLoad := int64(0)
	minLoad := int64(math.MaxInt64)

	for _, stats := range partStats {
		avgLoad += float64(stats.ComputeLoad)
		if stats.ComputeLoad > maxLoad {
			maxLoad = stats.ComputeLoad
		}
		if stats.ComputeLoad < minLoad {
			minLoad = stats.ComputeLoad
		}
	}
	avgLoad /= float64(nparts)

	imbalance := float64(maxLoad)/avgLoad - 1.0

	// Report statistics
	log.Printf("Partition Analysis:")
	log.Printf("  Objective value: %d", objval)
	log.Printf("  Cut edges: %d", cutEdges)
	log.Printf("  Communication volume: %d", commVolume)
	log.Printf("  Load imbalance: %.2f%%", imbalance*100)
	log.Printf("  Load range: [%d, %d], avg: %.1f", minLoad, maxLoad, avgLoad)

	// Report per-partition details
	log.Printf("\nPer-partition statistics:")
	for _, stats := range partStats {
		log.Printf("  Partition %d:", stats.ID)
		log.Printf("    Elements: %d", stats.NumElements)
		log.Printf("    Compute load: %d", stats.ComputeLoad)
		log.Printf("    Element types: %v", stats.ElementTypes)
		log.Printf("    Neighbors: %d", len(stats.NumNeighbors))
	}

	// Report interface statistics
	log.Printf("\nInterface statistics:")
	for pair, faces := range interfaceFaces {
		log.Printf("  Partition %d <-> %d: %d faces",
			pair[0], pair[1], len(faces))
	}
}

// PartitionStats holds statistics for a single partition
type PartitionStats struct {
	ID           int
	NumElements  int
	ComputeLoad  int64
	ElementTypes map[ElementType]int
	NumNeighbors map[int]int // neighbor partition -> shared faces
}

// Helper functions for partition analysis

// GetPartitionBoundaryFaces returns all faces on partition boundaries
func (mp *MeshPartitioner) GetPartitionBoundaryFaces() map[int][]int {
	boundaryFaces := make(map[int][]int) // partition -> face IDs

	for elem := 0; elem < mp.mesh.NumElements; elem++ {
		elemPart := mp.mesh.EToP[elem]

		for faceIdx, neighbor := range mp.mesh.EToE[elem] {
			isBoundary := false

			if neighbor < 0 {
				// Physical boundary
				isBoundary = true
			} else if mp.mesh.EToP[neighbor] != elemPart {
				// Partition boundary
				isBoundary = true
			}

			if isBoundary {
				faceID := mp.mesh.EToF[elem][faceIdx]
				boundaryFaces[elemPart] = append(boundaryFaces[elemPart], faceID)
			}
		}
	}

	return boundaryFaces
}

// GetPartitionElements returns all elements in a given partition
func (mp *MeshPartitioner) GetPartitionElements(partID int) []int {
	elements := []int{}
	for elem := 0; elem < mp.mesh.NumElements; elem++ {
		if mp.mesh.EToP[elem] == partID {
			elements = append(elements, elem)
		}
	}
	return elements
}

// ExportPartitionedMesh writes the partitioned mesh with partition information
func (mp *MeshPartitioner) ExportPartitionedMesh(filename string) error {
	// Implementation depends on desired output format
	// Could write VTK, Gmsh, or custom format with partition data
	return fmt.Errorf("not implemented")
}

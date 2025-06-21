package utils

import (
	"fmt"
	"strings"
	"unsafe"

	"github.com/notargets/gocca"
	"gonum.org/v1/gonum/mat"
)

// DataType represents the precision of numerical data
type DataType int

const (
	// Start from 1 to allow 0 to mean "unspecified"
	Float32 DataType = iota + 1
	Float64
	Int32
	Int64
)

// AlignmentType specifies memory alignment requirements
type AlignmentType int

const (
	NoAlignment    AlignmentType = 1
	CacheLineAlign AlignmentType = 64
	WarpAlign      AlignmentType = 128
	PageAlign      AlignmentType = 4096
)

// ArraySpec defines user requirements for array allocation
type ArraySpec struct {
	Name      string
	Size      int64         // Total size in bytes
	Alignment AlignmentType // Alignment requirement
	DataType  DataType      // Optional: specify Float32, Float64, Int32, or Int64
}

// arrayMetadata tracks information about allocated arrays
type arrayMetadata struct {
	spec             ArraySpec
	dataType         DataType // Float32, Float64, Int32, or Int64
	elementsPerValue int      // e.g., 1 for scalars, N for vectors
}

// KernelProgram manages code generation and execution for partition-parallel kernels
type KernelProgram struct {
	// Partition configuration
	NumPartitions int   // Number of partitions
	K             []int // Variable elements per partition

	// Type configuration
	FloatType DataType // Float32 or Float64
	IntType   DataType // Int32 or Int64

	// Static data to embed
	StaticMatrices map[string]mat.Matrix // Changed from Matrix to mat.Matrix

	// Array tracking for macro generation
	allocatedArrays []string
	arrayMetadata   map[string]arrayMetadata // Track array information

	// Generated code
	kernelPreamble string

	// Runtime resources
	device       *gocca.OCCADevice
	kernels      map[string]*gocca.OCCAKernel
	pooledMemory map[string]*gocca.OCCAMemory
}

// Config holds configuration for creating a KernelProgram
type Config struct {
	NumPartitions int   // Number of partitions
	K             []int // Elements per partition (variable)
	FloatType     DataType
	IntType       DataType
}

// NewKernelProgram creates a new KernelProgram with the given configuration
func NewKernelProgram(device *gocca.OCCADevice, cfg Config) *KernelProgram {
	// Validate device
	if device == nil {
		panic("device cannot be nil")
	}

	if !device.IsInitialized() {
		panic("device is not initialized")
	}

	// Validate partition configuration
	if len(cfg.K) == 0 {
		panic("K array cannot be empty")
	}

	if cfg.NumPartitions == 0 {
		cfg.NumPartitions = len(cfg.K)
	} else if cfg.NumPartitions != len(cfg.K) {
		panic("NumPartitions must match length of K array")
	}

	// Set type defaults only if both are zero
	if cfg.FloatType == 0 && cfg.IntType == 0 {
		cfg.FloatType = Float64
		cfg.IntType = Int64
	}

	kp := &KernelProgram{
		NumPartitions:   cfg.NumPartitions,
		K:               make([]int, len(cfg.K)),
		FloatType:       cfg.FloatType,
		IntType:         cfg.IntType,
		StaticMatrices:  make(map[string]mat.Matrix), // Changed from Matrix to mat.Matrix
		allocatedArrays: []string{},
		arrayMetadata:   make(map[string]arrayMetadata),
		device:          device,
		kernels:         make(map[string]*gocca.OCCAKernel),
		pooledMemory:    make(map[string]*gocca.OCCAMemory),
	}

	// Copy K values
	copy(kp.K, cfg.K)

	return kp
}

// AddStaticMatrix adds a matrix that will be compiled into kernels as static data
func (kp *KernelProgram) AddStaticMatrix(name string, matrix mat.Matrix) { // Changed from Matrix to mat.Matrix
	kp.StaticMatrices[name] = matrix
}

// AllocateArrays allocates device memory with user-specified requirements
func (kp *KernelProgram) AllocateArrays(specs []ArraySpec) error {
	for _, spec := range specs {
		// Step 1: Calculate aligned offsets first
		offsets, totalSize := kp.calculateAlignedOffsetsAndSize(spec)

		// Step 2: Allocate global pool with exact calculated size
		memory := kp.device.Malloc(totalSize, nil, nil)
		if memory == nil {
			return fmt.Errorf("failed to allocate %s", spec.Name)
		}
		kp.pooledMemory[spec.Name+"_global"] = memory

		// Step 3: Store offsets array - convert to appropriate int type
		offsetSize := int64(len(offsets)) * int64(kp.GetIntSize())
		var offsetMem *gocca.OCCAMemory

		if kp.GetIntSize() == 4 {
			// Convert int64 offsets to int32 before storing
			offsets32 := make([]int32, len(offsets))
			for i, v := range offsets {
				offsets32[i] = int32(v)
			}
			offsetMem = kp.device.Malloc(offsetSize, unsafe.Pointer(&offsets32[0]), nil)
		} else {
			// Store int64 offsets directly
			offsetMem = kp.device.Malloc(offsetSize, unsafe.Pointer(&offsets[0]), nil)
		}

		if offsetMem == nil {
			return fmt.Errorf("failed to allocate offsets for %s", spec.Name)
		}
		kp.pooledMemory[spec.Name+"_offsets"] = offsetMem

		// Track array name for macro generation
		kp.allocatedArrays = append(kp.allocatedArrays, spec.Name)

		// Store metadata for later retrieval
		// Use specified data type, or default to FloatType
		dataType := spec.DataType
		if dataType == 0 {
			dataType = kp.FloatType
		}

		// Calculate element size based on data type
		var elementSize int
		switch dataType {
		case Float32, Int32:
			elementSize = 4
		case Float64, Int64:
			elementSize = 8
		default:
			elementSize = kp.GetFloatSize()
		}

		elementsPerValue := int(spec.Size) / (kp.getTotalElements() * elementSize)

		kp.arrayMetadata[spec.Name] = arrayMetadata{
			spec:             spec,
			dataType:         dataType,
			elementsPerValue: elementsPerValue,
		}
	}

	// Allocate K array - same fix needed here
	kSize := int64(len(kp.K)) * int64(kp.GetIntSize())
	var kMem *gocca.OCCAMemory

	if kp.GetIntSize() == 4 {
		// Convert int K values to int32
		k32 := make([]int32, len(kp.K))
		for i, v := range kp.K {
			k32[i] = int32(v)
		}
		kMem = kp.device.Malloc(kSize, unsafe.Pointer(&k32[0]), nil)
	} else {
		// Convert int K values to int64
		k64 := make([]int64, len(kp.K))
		for i, v := range kp.K {
			k64[i] = int64(v)
		}
		kMem = kp.device.Malloc(kSize, unsafe.Pointer(&k64[0]), nil)
	}

	if kMem == nil {
		return fmt.Errorf("failed to allocate K array")
	}
	kp.pooledMemory["K"] = kMem

	return nil
}

// calculateAlignedOffsetsAndSize calculates partition offsets and total size needed
func (kp *KernelProgram) calculateAlignedOffsetsAndSize(spec ArraySpec) ([]int64, int64) {
	offsets := make([]int64, kp.NumPartitions+1)

	// Calculate bytes per element based on spec size and total elements
	totalElements := kp.getTotalElements()
	bytesPerElement := spec.Size / int64(totalElements)

	alignment := int64(spec.Alignment)
	currentByteOffset := int64(0)

	for i := 0; i < kp.NumPartitions; i++ {
		// Align current offset
		if currentByteOffset%alignment != 0 {
			currentByteOffset = ((currentByteOffset + alignment - 1) / alignment) * alignment
		}

		// Store element offset (not byte offset)
		offsets[i] = currentByteOffset / bytesPerElement

		// Advance by partition data size
		partitionBytes := int64(kp.K[i]) * bytesPerElement
		currentByteOffset += partitionBytes
	}

	// Final offset for bounds checking
	if currentByteOffset%alignment != 0 {
		currentByteOffset = ((currentByteOffset + alignment - 1) / alignment) * alignment
	}
	offsets[kp.NumPartitions] = currentByteOffset / bytesPerElement

	// Total size is the final byte offset
	totalSize := currentByteOffset

	return offsets, totalSize
}

// getTotalElements returns sum of all K values
func (kp *KernelProgram) getTotalElements() int {
	total := 0
	for _, k := range kp.K {
		total += k
	}
	return total
}

// GeneratePreamble generates the kernel preamble with static data and utilities
func (kp *KernelProgram) GeneratePreamble() string {
	var sb strings.Builder

	// 1. Type definitions and constants
	sb.WriteString(kp.generateTypeDefinitions())

	// 2. Static matrix declarations
	sb.WriteString(kp.generateStaticMatrices())

	// 3. Partition access macros
	sb.WriteString(kp.generatePartitionMacros())

	// 4. Matrix operation macros
	sb.WriteString(kp.generateMatrixMacros())

	kp.kernelPreamble = sb.String()
	return kp.kernelPreamble
}

// generateTypeDefinitions creates type definitions based on precision settings
func (kp *KernelProgram) generateTypeDefinitions() string {
	var sb strings.Builder

	// Type definitions
	floatTypeStr := "double"
	floatSuffix := ""
	if kp.FloatType == Float32 {
		floatTypeStr = "float"
		floatSuffix = "f"
	}

	intTypeStr := "long"
	if kp.IntType == Int32 {
		intTypeStr = "int"
	}

	sb.WriteString(fmt.Sprintf("typedef %s real_t;\n", floatTypeStr))
	sb.WriteString(fmt.Sprintf("typedef %s int_t;\n", intTypeStr))
	sb.WriteString(fmt.Sprintf("#define REAL_ZERO 0.0%s\n", floatSuffix))
	sb.WriteString(fmt.Sprintf("#define REAL_ONE 1.0%s\n", floatSuffix))
	sb.WriteString("\n")

	// Constants - these would typically be set based on the specific kernel needs
	sb.WriteString(fmt.Sprintf("#define NPART %d\n", kp.NumPartitions))
	sb.WriteString("\n")

	return sb.String()
}

// generateStaticMatrices converts matrices to static array initializations
func (kp *KernelProgram) generateStaticMatrices() string {
	var sb strings.Builder

	if len(kp.StaticMatrices) > 0 {
		sb.WriteString("// Static matrices\n")
		for name, matrix := range kp.StaticMatrices {
			sb.WriteString(kp.formatStaticMatrix(name, matrix))
		}
		sb.WriteString("\n")
	}

	return sb.String()
}

// formatStaticMatrix formats a single matrix as a static C array
func (kp *KernelProgram) formatStaticMatrix(name string, m mat.Matrix) string { // Changed from Matrix to mat.Matrix
	rows, cols := m.Dims()
	var sb strings.Builder

	typeStr := "double"
	if kp.FloatType == Float32 {
		typeStr = "float"
	}

	sb.WriteString(fmt.Sprintf("const %s %s[%d][%d] = {\n",
		typeStr, name, rows, cols))

	for i := 0; i < rows; i++ {
		sb.WriteString("    {")
		for j := 0; j < cols; j++ {
			if j > 0 {
				sb.WriteString(", ")
			}
			val := m.At(i, j)
			if kp.FloatType == Float32 {
				sb.WriteString(fmt.Sprintf("%.7ef", val))
			} else {
				sb.WriteString(fmt.Sprintf("%.15e", val))
			}
		}
		sb.WriteString("}")
		if i < rows-1 {
			sb.WriteString(",")
		}
		sb.WriteString("\n")
	}
	sb.WriteString("};\n\n")

	return sb.String()
}

// generatePartitionMacros creates macros for partition data access
func (kp *KernelProgram) generatePartitionMacros() string {
	var sb strings.Builder

	sb.WriteString("// Partition access macros\n")

	// Generate macro for each allocated array
	for _, arrayName := range kp.allocatedArrays {
		sb.WriteString(fmt.Sprintf("#define %s_PART(part) (%s_global + %s_offsets[part])\n",
			arrayName, arrayName, arrayName))
	}

	if len(kp.allocatedArrays) > 0 {
		sb.WriteString("\n")
	}

	return sb.String()
}

// generateMatrixMacros creates vectorizable matrix multiplication macros
func (kp *KernelProgram) generateMatrixMacros() string {
	var sb strings.Builder

	sb.WriteString("// Vectorizable matrix multiplication macros\n")

	// Generate macro for each static matrix
	for name, matrix := range kp.StaticMatrices {
		rows, cols := matrix.Dims()

		// Determine if it's a square matrix
		if rows == cols {
			// Square matrix macro (like differentiation matrices)
			sb.WriteString(fmt.Sprintf("#define MATMUL_%s(IN, OUT, K_VAL, NP) do { \\\n", name))
			sb.WriteString("    for (int i = 0; i < (NP); ++i) { \\\n")
			sb.WriteString("        for (int elem = 0; elem < (K_VAL); ++elem) { \\\n")
			sb.WriteString("            real_t sum = REAL_ZERO; \\\n")
			sb.WriteString("            for (int j = 0; j < (NP); ++j) { \\\n")
			sb.WriteString(fmt.Sprintf("                sum += %s[i][j] * (IN)[elem * (NP) + j]; \\\n", name))
			sb.WriteString("            } \\\n")
			sb.WriteString("            (OUT)[elem * (NP) + i] = sum; \\\n")
			sb.WriteString("        } \\\n")
			sb.WriteString("    } \\\n")
			sb.WriteString("} while(0)\n\n")
		} else {
			// Rectangular matrix macro (need to know input/output dimensions)
			// For now, generate a comment explaining what's needed
			sb.WriteString(fmt.Sprintf("// Note: %s is a %dx%d rectangular matrix\n", name, rows, cols))
			sb.WriteString(fmt.Sprintf("// Macro generation requires knowledge of input dimensions\n\n"))
		}
	}

	return sb.String()
}

// BuildKernel compiles and registers a kernel with the program
func (kp *KernelProgram) BuildKernel(kernelSource, kernelName string) (*gocca.OCCAKernel, error) {
	// Generate preamble if not already done
	if kp.kernelPreamble == "" {
		kp.GeneratePreamble()
	}

	// Combine preamble with kernel source
	fullSource := kp.kernelPreamble + "\n" + kernelSource

	// Build kernel
	kernel, err := kp.device.BuildKernelFromString(fullSource, kernelName, nil)
	if err != nil {
		return nil, fmt.Errorf("failed to build kernel %s: %w", kernelName, err)
	}

	// Register if build succeeded
	if kernel != nil {
		kp.kernels[kernelName] = kernel
		return kernel, nil
	}

	return nil, fmt.Errorf("kernel build returned nil for %s", kernelName)
}

// RunKernel executes a registered kernel with the given arguments
func (kp *KernelProgram) RunKernel(name string, args ...interface{}) error {
	kernel, exists := kp.kernels[name]
	if !exists {
		return fmt.Errorf("kernel %s not found", name)
	}

	// Configure for partition-parallel execution
	outerDims := gocca.OCCADim{
		X: uint64(kp.NumPartitions),
		Y: 1,
		Z: 1,
	}

	// Inner dimensions would be set based on the specific kernel requirements
	// For now, use a default
	innerDims := gocca.OCCADim{
		X: 32, // Default work group size
		Y: 1,
		Z: 1,
	}

	kernel.SetRunDims(outerDims, innerDims)

	// Expand args to include renamed arrays
	expandedArgs := kp.expandKernelArgs(args)

	return kernel.RunWithArgs(expandedArgs...)
}

// expandKernelArgs transforms user array names to kernel parameter names
func (kp *KernelProgram) expandKernelArgs(args []interface{}) []interface{} {
	expanded := []interface{}{}

	// Always pass K array first
	expanded = append(expanded, kp.pooledMemory["K"])

	// Process remaining arguments
	for _, arg := range args {
		switch v := arg.(type) {
		case string:
			// Check if this is an allocated array name
			if kp.isAllocatedArray(v) {
				// Add both _global and _offsets
				expanded = append(expanded, kp.pooledMemory[v+"_global"])
				expanded = append(expanded, kp.pooledMemory[v+"_offsets"])
			} else {
				// Not an array name, pass as-is
				expanded = append(expanded, arg)
			}
		default:
			// Pass through non-string arguments
			expanded = append(expanded, arg)
		}
	}

	return expanded
}

// isAllocatedArray checks if a name corresponds to an allocated array
func (kp *KernelProgram) isAllocatedArray(name string) bool {
	for _, arrayName := range kp.allocatedArrays {
		if arrayName == name {
			return true
		}
	}
	return false
}

// GetMemory returns a device memory handle by name
func (kp *KernelProgram) GetMemory(name string) *gocca.OCCAMemory {
	// First check if it's a direct memory name
	if mem, exists := kp.pooledMemory[name]; exists {
		return mem
	}

	// Check if it's an allocated array (return the _global version)
	if kp.isAllocatedArray(name) {
		return kp.pooledMemory[name+"_global"]
	}

	return nil
}

// GetFloatSize returns the size in bytes of the float type
func (kp *KernelProgram) GetFloatSize() int {
	if kp.FloatType == Float32 {
		return 4
	}
	return 8
}

// GetIntSize returns the size in bytes of the int type
func (kp *KernelProgram) GetIntSize() int {
	if kp.IntType == Int32 {
		return 4
	}
	return 8
}

// Free releases all allocated resources
func (kp *KernelProgram) Free() {
	// Free kernels
	for _, kernel := range kp.kernels {
		if kernel != nil {
			kernel.Free()
		}
	}

	// Free memory
	for _, mem := range kp.pooledMemory {
		if mem != nil {
			mem.Free()
		}
	}
}

// GetArrayType returns the data type of an allocated array
func (kp *KernelProgram) GetArrayType(name string) (DataType, error) {
	metadata, exists := kp.arrayMetadata[name]
	if !exists {
		return 0, fmt.Errorf("array %s not found", name)
	}
	return metadata.dataType, nil
}

// GetArrayLogicalSize returns the number of logical elements in an array
func (kp *KernelProgram) GetArrayLogicalSize(name string) (int, error) {
	metadata, exists := kp.arrayMetadata[name]
	if !exists {
		return 0, fmt.Errorf("array %s not found", name)
	}

	// Total elements = total elements across all partitions * elements per value
	return kp.getTotalElements() * metadata.elementsPerValue, nil
}

// CopyArrayToHost copies array data from device to host, removing alignment padding
func CopyArrayToHost[T any](kp *KernelProgram, name string) ([]T, error) {
	// Check if array exists
	metadata, exists := kp.arrayMetadata[name]
	if !exists {
		return nil, fmt.Errorf("array %s not found", name)
	}

	// Verify type matches
	var sample T
	requestedType := getDataTypeFromSample(sample)
	if requestedType != metadata.dataType {
		return nil, fmt.Errorf("type mismatch: array is %v, requested %v",
			metadata.dataType, requestedType)
	}

	// Get memory and offsets
	memory := kp.GetMemory(name)
	if memory == nil {
		return nil, fmt.Errorf("memory for %s not found", name)
	}

	offsetsMem := kp.pooledMemory[name+"_offsets"]
	if offsetsMem == nil {
		return nil, fmt.Errorf("offsets for %s not found", name)
	}

	// Read offsets to determine actual data locations
	numOffsets := kp.NumPartitions + 1
	var offsets []int64

	if kp.GetIntSize() == 4 {
		offsets32 := make([]int32, numOffsets)
		offsetsMem.CopyTo(unsafe.Pointer(&offsets32[0]), int64(numOffsets*4))
		offsets = make([]int64, numOffsets)
		for i, v := range offsets32 {
			offsets[i] = int64(v)
		}
	} else {
		offsets = make([]int64, numOffsets)
		offsetsMem.CopyTo(unsafe.Pointer(&offsets[0]), int64(numOffsets*8))
	}

	// Calculate total logical elements
	totalLogicalElements := kp.getTotalElements() * metadata.elementsPerValue

	// Allocate result array
	result := make([]T, totalLogicalElements)

	// Determine element size
	elementSize := int64(unsafe.Sizeof(sample))

	// Read all data including padding
	totalSize := offsets[kp.NumPartitions] * elementSize
	tempData := make([]T, offsets[kp.NumPartitions])
	memory.CopyTo(unsafe.Pointer(&tempData[0]), totalSize)

	// Copy data partition by partition, skipping alignment padding
	destIdx := 0
	for part := 0; part < kp.NumPartitions; part++ {
		startIdx := offsets[part]
		partitionElements := kp.K[part] * metadata.elementsPerValue

		for i := 0; i < partitionElements; i++ {
			result[destIdx] = tempData[startIdx+int64(i)]
			destIdx++
		}
	}

	return result, nil
}

// CopyPartitionToHost copies a specific partition's data from device to host
func CopyPartitionToHost[T any](kp *KernelProgram, name string, partitionID int) ([]T, error) {
	if partitionID < 0 || partitionID >= kp.NumPartitions {
		return nil, fmt.Errorf("invalid partition ID %d", partitionID)
	}

	// Check if array exists
	metadata, exists := kp.arrayMetadata[name]
	if !exists {
		return nil, fmt.Errorf("array %s not found", name)
	}

	// Verify type matches
	var sample T
	requestedType := getDataTypeFromSample(sample)
	if requestedType != metadata.dataType {
		return nil, fmt.Errorf("type mismatch: array is %v, requested %v",
			metadata.dataType, requestedType)
	}

	// Calculate partition size
	partitionElements := kp.K[partitionID] * metadata.elementsPerValue

	// Handle empty partition case early
	if partitionElements == 0 {
		return make([]T, 0), nil
	}

	// Get memory and offsets
	globalMem := kp.pooledMemory[name+"_global"]
	offsetsMem := kp.pooledMemory[name+"_offsets"]

	if globalMem == nil || offsetsMem == nil {
		return nil, fmt.Errorf("memory not allocated for array %s", name)
	}

	// Get offsets - fixed to handle int32 correctly
	numOffsets := kp.NumPartitions + 1
	offsets := make([]int64, numOffsets)

	if kp.GetIntSize() == 4 {
		// For int32, copy to int32 array first, then convert
		offsets32 := make([]int32, numOffsets)
		offsetSize := int64(numOffsets * 4)
		offsetsMem.CopyTo(unsafe.Pointer(&offsets32[0]), offsetSize)
		for i, v := range offsets32 {
			offsets[i] = int64(v)
		}
	} else {
		// For int64, copy directly
		offsetSize := int64(numOffsets * 8)
		offsetsMem.CopyTo(unsafe.Pointer(&offsets[0]), offsetSize)
	}

	// Allocate result array
	result := make([]T, partitionElements)

	// Determine element size
	elementSize := int64(unsafe.Sizeof(sample))

	// Read all data including padding
	totalSize := offsets[kp.NumPartitions] * elementSize
	tempData := make([]T, offsets[kp.NumPartitions])
	globalMem.CopyTo(unsafe.Pointer(&tempData[0]), totalSize)

	// Copy just this partition's data
	startIdx := offsets[partitionID]
	for i := 0; i < partitionElements; i++ {
		result[i] = tempData[startIdx+int64(i)]
	}

	return result, nil
}

// Helper function to determine DataType from a sample value
func getDataTypeFromSample(sample interface{}) DataType {
	switch sample.(type) {
	case float32:
		return Float32
	case float64:
		return Float64
	case int32:
		return Int32
	case int64:
		return Int64
	default:
		return 0
	}
}

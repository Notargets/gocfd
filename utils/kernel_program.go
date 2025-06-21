package utils

import (
	"fmt"
	"strings"
	"unsafe"

	"github.com/notargets/gocca"
)

// DataType represents the precision of numerical data
type DataType int

const (
	Float32 DataType = iota
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
	StaticMatrices map[string]Matrix

	// Array tracking for macro generation
	allocatedArrays []string

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
		StaticMatrices:  make(map[string]Matrix),
		allocatedArrays: []string{},
		device:          device,
		kernels:         make(map[string]*gocca.OCCAKernel),
		pooledMemory:    make(map[string]*gocca.OCCAMemory),
	}

	// Copy K values
	copy(kp.K, cfg.K)

	return kp
}

// AddStaticMatrix adds a matrix that will be compiled into kernels as static data
func (kp *KernelProgram) AddStaticMatrix(name string, matrix Matrix) {
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

		// Step 3: Store offsets array
		offsetSize := int64(len(offsets)) * int64(kp.GetIntSize())
		offsetMem := kp.device.Malloc(offsetSize, unsafe.Pointer(&offsets[0]), nil)
		if offsetMem == nil {
			return fmt.Errorf("failed to allocate offsets for %s", spec.Name)
		}
		kp.pooledMemory[spec.Name+"_offsets"] = offsetMem

		// Track array name for macro generation
		kp.allocatedArrays = append(kp.allocatedArrays, spec.Name)
	}

	// Allocate K array
	kSize := int64(len(kp.K)) * int64(kp.GetIntSize())
	kMem := kp.device.Malloc(kSize, unsafe.Pointer(&kp.K[0]), nil)
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
func (kp *KernelProgram) formatStaticMatrix(name string, m Matrix) string {
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
			sb.WriteString(fmt.Sprintf(`#define MATMUL_%s(IN, OUT, K_VAL, NP) do { \
    for (int i = 0; i < (NP); ++i) { \
        for (int elem = 0; elem < (K_VAL); ++elem) { \
            real_t sum = REAL_ZERO; \
            for (int j = 0; j < (NP); ++j) { \
                sum += %s[i][j] * (IN)[elem * (NP) + j]; \
            } \
            (OUT)[elem * (NP) + i] = sum; \
        } \
    } \
} while(0)

`, name, name))
		} else {
			// Rectangular matrix (like LIFT)
			sb.WriteString(fmt.Sprintf(`#define MATMUL_%s(IN, OUT, K_VAL, ROWS, IN_COLS) do { \
    for (int i = 0; i < (ROWS); ++i) { \
        for (int elem = 0; elem < (K_VAL); ++elem) { \
            real_t sum = REAL_ZERO; \
            for (int j = 0; j < (IN_COLS); ++j) { \
                sum += %s[i][j] * (IN)[elem * (IN_COLS) + j]; \
            } \
            (OUT)[elem * (ROWS) + i] = sum; \
        } \
    } \
} while(0)

`, name, name))
		}
	}

	return sb.String()
}

// BuildKernel compiles a kernel from source with the generated preamble
func (kp *KernelProgram) BuildKernel(kernelSource, kernelName string) (*gocca.OCCAKernel, error) {
	if kp.device == nil {
		return nil, fmt.Errorf("device is nil")
	}

	if !kp.device.IsInitialized() {
		return nil, fmt.Errorf("device is not initialized")
	}

	if kernelName == "" {
		return nil, fmt.Errorf("kernel name cannot be empty")
	}

	if kernelSource == "" {
		return nil, fmt.Errorf("kernel source cannot be empty")
	}

	// Ensure preamble is generated
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

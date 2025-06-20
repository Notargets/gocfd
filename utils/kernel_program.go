package utils

import (
	"fmt"
	"strings"

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

// KernelProgram manages code generation and execution for DG kernels
type KernelProgram struct {
	// Configuration
	Order       int      // Polynomial order (N)
	Np          int      // Number of nodes per element
	Nfp         int      // Number of face points
	NumElements int      // Number of elements in this partition (K)
	FloatType   DataType // Float32 or Float64 (default: Float64)
	IntType     DataType // Int32 or Int64 (default: Int64)

	// Static data to embed
	StaticMatrices map[string]Matrix // Dr, Ds, Dt, LIFT, etc.

	// Generated code
	kernelPreamble string // Generated static data and utilities

	// Runtime resources
	device  *gocca.OCCADevice
	kernels map[string]*gocca.OCCAKernel
	memory  map[string]*gocca.OCCAMemory
}

// Config holds configuration for creating a KernelProgram
type Config struct {
	Order       int
	NumElements int
	FloatType   DataType
	IntType     DataType
}

// New creates a new KernelProgram with the given configuration
func NewKernelProgram(device *gocca.OCCADevice, cfg Config) *KernelProgram {
	// Set defaults only if both are zero
	if cfg.FloatType == 0 && cfg.IntType == 0 {
		cfg.FloatType = Float64
		cfg.IntType = Int64
	}

	// Calculate derived values
	np := (cfg.Order + 1) * (cfg.Order + 2) * (cfg.Order + 3) / 6
	nfp := (cfg.Order + 1) * (cfg.Order + 2) / 2

	kp := &KernelProgram{
		Order:          cfg.Order,
		Np:             np,
		Nfp:            nfp,
		NumElements:    cfg.NumElements,
		FloatType:      cfg.FloatType,
		IntType:        cfg.IntType,
		StaticMatrices: make(map[string]Matrix),
		device:         device,
		kernels:        make(map[string]*gocca.OCCAKernel),
		memory:         make(map[string]*gocca.OCCAMemory),
	}

	return kp
}

// AddStaticMatrix adds a matrix that will be compiled into kernels as static data
func (kp *KernelProgram) AddStaticMatrix(name string, matrix Matrix) {
	kp.StaticMatrices[name] = matrix
}

// GenerateKernelMain generates the kernel preamble with static data and utilities
func (kp *KernelProgram) GenerateKernelMain() string {
	var sb strings.Builder

	// 1. Type definitions and constants
	sb.WriteString(kp.generateTypeDefinitions())

	// 2. Static matrix declarations
	sb.WriteString(kp.generateStaticMatrices())

	// 3. Utility functions
	sb.WriteString(kp.generateUtilityFunctions())

	kp.kernelPreamble = sb.String()
	return kp.kernelPreamble
}

// generateTypeDefinitions creates type definitions based on precision settings
func (kp *KernelProgram) generateTypeDefinitions() string {
	var sb strings.Builder

	// Constants
	sb.WriteString(fmt.Sprintf("#define ORDER %d\n", kp.Order))
	sb.WriteString(fmt.Sprintf("#define NP %d\n", kp.Np))
	sb.WriteString(fmt.Sprintf("#define NFP %d\n", kp.Nfp))
	sb.WriteString("\n")

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

	// Face type constants
	sb.WriteString("// Face types\n")
	sb.WriteString("#define FACE_INTERIOR 0\n")
	sb.WriteString("#define FACE_BOUNDARY 1\n")
	sb.WriteString("\n")

	// Boundary condition types
	sb.WriteString("// Boundary condition types\n")
	sb.WriteString("#define BC_WALL 1\n")
	sb.WriteString("#define BC_INFLOW 2\n")
	sb.WriteString("#define BC_OUTFLOW 3\n")
	sb.WriteString("\n")

	return sb.String()
}

// generateStaticMatrices converts matrices to static array initializations
func (kp *KernelProgram) generateStaticMatrices() string {
	var sb strings.Builder

	sb.WriteString("// Static matrices\n")

	// Generate each matrix
	for name, matrix := range kp.StaticMatrices {
		sb.WriteString(kp.formatStaticMatrix(name, matrix))
	}

	sb.WriteString("\n")
	return sb.String()
}

// formatStaticMatrix formats a single matrix as a static C array
func (kp *KernelProgram) formatStaticMatrix(name string, m Matrix) string {
	rows, cols := m.Dims()
	var sb strings.Builder

	// Use appropriate type based on FloatType
	typeStr := "double"
	if kp.FloatType == Float32 {
		typeStr = "float"
	}

	// Use const instead of __constant__ for OCCA
	sb.WriteString(fmt.Sprintf("const %s %s[%d][%d] = {\n",
		typeStr, name, rows, cols))

	for i := 0; i < rows; i++ {
		sb.WriteString("    {")
		for j := 0; j < cols; j++ {
			if j > 0 {
				sb.WriteString(", ")
			}
			// Use appropriate precision
			if kp.FloatType == Float32 {
				sb.WriteString(fmt.Sprintf("%.7ef", m.At(i, j)))
			} else {
				sb.WriteString(fmt.Sprintf("%.15e", m.At(i, j)))
			}
		}
		sb.WriteString("},\n")
	}
	sb.WriteString("};\n\n")

	return sb.String()
}

// generateUtilityFunctions creates utility functions for matrix operations
func (kp *KernelProgram) generateUtilityFunctions() string {
	var sb strings.Builder

	sb.WriteString("// Utility functions\n\n")

	// Matrix multiplication functions for each static matrix
	for name := range kp.StaticMatrices {
		sb.WriteString(kp.generateMatMulFunction(name))
	}

	// Generic utilities
	sb.WriteString(kp.generateGenericUtilities())

	return sb.String()
}

// generateMatMulFunction creates a matrix multiplication function for a specific static matrix
func (kp *KernelProgram) generateMatMulFunction(matrixName string) string {
	return fmt.Sprintf(`// Matrix multiplication using static %s
inline void matMul_%s_Large(const real_t* U, real_t* result, int K) {
    for (int i = 0; i < NP; ++i) {
        for (int k = 0; k < K; ++k) {
            real_t sum = REAL_ZERO;
            #pragma unroll
            for (int j = 0; j < NP; ++j) {
                sum += %s[i][j] * U[j*K + k];
            }
            result[i*K + k] = sum;
        }
    }
}

`, matrixName, matrixName, matrixName)
}

// generateGenericUtilities creates general-purpose utility functions
func (kp *KernelProgram) generateGenericUtilities() string {
	return `// Generic matrix-vector multiplication
inline void matVecMul(const real_t* mat, const real_t* vec, real_t* result, int M, int N) {
    for (int i = 0; i < M; ++i) {
        real_t sum = REAL_ZERO;
        #pragma unroll
        for (int j = 0; j < N; ++j) {
            sum += mat[i*N + j] * vec[j];
        }
        result[i] = sum;
    }
}

// Apply boundary condition
inline real_t applyBC(int bcType, real_t M, real_t bcData,
                     real_t nx, real_t ny, real_t nz) {
    switch(bcType) {
        case BC_WALL:
            return -M;  // Reflection
        case BC_INFLOW:
            return bcData;  // Prescribed value
        case BC_OUTFLOW:
            return M;  // Zero gradient
        default:
            return M;
    }
}

// Simple Riemann flux (Rusanov/Lax-Friedrichs)
inline real_t riemannFlux(real_t M, real_t P, 
                         real_t nx, real_t ny, real_t nz) {
    real_t jump = P - M;
    real_t wavespeed = REAL_ONE; // Simplified - would compute actual wave speed
    return (real_t)0.5 * ((M + P) - wavespeed * jump);
}

`
}

// AllocateKernelMemory allocates device memory for runtime data
func (kp *KernelProgram) AllocateKernelMemory() error {
	// Get size based on precision
	floatSize := 8 // Float64
	if kp.FloatType == Float32 {
		floatSize = 4
	}

	intSize := 8 // Int64
	if kp.IntType == Int32 {
		intSize = 4
	}

	nodeCount := kp.Np * kp.NumElements
	faceCount := 4 * kp.Nfp * kp.NumElements

	// Solution arrays
	kp.memory["U"] = kp.device.Malloc(
		int64(nodeCount*floatSize),
		nil,
		nil,
	)

	// RHS array
	kp.memory["RHS"] = kp.device.Malloc(
		int64(nodeCount*floatSize),
		nil,
		nil,
	)

	// Geometric factors (stored separately for vectorization)
	geoFactors := []string{"rx", "ry", "rz", "sx", "sy", "sz", "tx", "ty", "tz"}
	for _, factor := range geoFactors {
		kp.memory[factor] = kp.device.Malloc(
			int64(nodeCount*floatSize),
			nil,
			nil,
		)
	}

	// Face data
	kp.memory["faceM"] = kp.device.Malloc(
		int64(faceCount*floatSize),
		nil,
		nil,
	)

	kp.memory["faceP"] = kp.device.Malloc(
		int64(faceCount*floatSize),
		nil,
		nil,
	)

	// Face types (integer)
	kp.memory["faceTypes"] = kp.device.Malloc(
		int64(faceCount*intSize),
		nil,
		nil,
	)

	// Face normals
	faceNormals := []string{"nx", "ny", "nz"}
	for _, normal := range faceNormals {
		kp.memory[normal] = kp.device.Malloc(
			int64(faceCount*floatSize),
			nil,
			nil,
		)
	}

	// Face scaling
	kp.memory["Fscale"] = kp.device.Malloc(
		int64(faceCount*floatSize),
		nil,
		nil,
	)

	// Boundary condition data
	kp.memory["bcData"] = kp.device.Malloc(
		int64(faceCount*floatSize),
		nil,
		nil,
	)

	return nil
}

// RegisterKernel adds a compiled kernel to the program
func (kp *KernelProgram) RegisterKernel(name string, kernel *gocca.OCCAKernel) {
	kp.kernels[name] = kernel
}

// BuildKernel compiles a kernel from source with the generated preamble
func (kp *KernelProgram) BuildKernel(kernelSource, kernelName string) (*gocca.OCCAKernel, error) {
	// Combine preamble with kernel source
	fullSource := kp.kernelPreamble + "\n" + kernelSource

	// Build kernel
	kernel, err := kp.device.BuildKernelFromString(fullSource, kernelName, nil)
	if err != nil {
		return nil, fmt.Errorf("failed to build kernel %s: %w", kernelName, err)
	}

	// Register it
	kp.RegisterKernel(kernelName, kernel)

	return kernel, nil
}

// RunKernel executes a registered kernel with the given arguments
func (kp *KernelProgram) RunKernel(name string, args ...interface{}) error {
	kernel, exists := kp.kernels[name]
	if !exists {
		return fmt.Errorf("kernel %s not found", name)
	}

	// Single kernel execution for entire partition
	// No thread dimensions needed - kernel internally vectorizes
	kernel.SetRunDims(
		gocca.OCCADim{X: 1, Y: 1, Z: 1}, // Single execution
		gocca.OCCADim{X: 1, Y: 1, Z: 1},
	)

	return kernel.RunWithArgs(args...)
}

// GetMemory returns a device memory handle by name
func (kp *KernelProgram) GetMemory(name string) *gocca.OCCAMemory {
	return kp.memory[name]
}

// Free releases all allocated resources
func (kp *KernelProgram) Free() {
	// Free kernels
	for _, kernel := range kp.kernels {
		kernel.Free()
	}

	// Free memory
	for _, mem := range kp.memory {
		mem.Free()
	}
}

// GetKernelPreamble returns the generated preamble (useful for debugging)
func (kp *KernelProgram) GetKernelPreamble() string {
	return kp.kernelPreamble
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

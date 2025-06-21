package utils

import (
	"fmt"
	"math"
	"strings"
	"testing"
	"unsafe"

	"github.com/notargets/gocca"
)

// Test 1: Creation and Configuration (Fundamental)
func TestKernelProgram_Creation(t *testing.T) {
	// Test basic object creation with different orders
	for _, order := range []int{1, 2, 3, 4} {
		t.Run(fmt.Sprintf("Order_%d", order), func(t *testing.T) {
			device := createTestDevice(t)
			defer device.Free()

			kp := NewKernelProgram(device, Config{
				Order:           order,
				NumPartitions:   1,
				ElementsPerPart: 10,
				// Don't set types - let them default
			})
			defer kp.Free()

			// Validate Np calculation
			expectedNp := (order + 1) * (order + 2) * (order + 3) / 6
			if kp.Np != expectedNp {
				t.Errorf("Order %d: Expected Np=%d, got %d",
					order, expectedNp, kp.Np)
			}

			// Validate Nfp calculation
			expectedNfp := (order + 1) * (order + 2) / 2
			if kp.Nfp != expectedNfp {
				t.Errorf("Order %d: Expected Nfp=%d, got %d",
					order, expectedNfp, kp.Nfp)
			}

			// Verify default types
			if kp.FloatType != Float64 {
				t.Errorf("Expected default FloatType=Float64, got %v", kp.FloatType)
			}
			if kp.IntType != Int64 {
				t.Errorf("Expected default IntType=Int64, got %v", kp.IntType)
			}

			// Verify default partition count
			if kp.NumPartitions != 1 {
				t.Errorf("Expected default NumPartitions=1, got %d", kp.NumPartitions)
			}
		})
	}
}

// Test partition configurations
func TestKernelProgram_PartitionConfigurations(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	testCases := []struct {
		name            string
		numPartitions   int
		elementsPerPart int
	}{
		{"SinglePartition", 0, 100}, // 0 should default to 1
		{"TwoPartitions", 2, 50},
		{"EightPartitions", 8, 25},
		{"ManyPartitions", 64, 10},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			kp := NewKernelProgram(device, Config{
				Order:           3,
				NumPartitions:   tc.numPartitions,
				ElementsPerPart: tc.elementsPerPart,
			})
			defer kp.Free()

			expectedPartitions := tc.numPartitions
			if expectedPartitions == 0 {
				expectedPartitions = 1 // Default
			}

			if kp.NumPartitions != expectedPartitions {
				t.Errorf("Expected NumPartitions=%d, got %d",
					expectedPartitions, kp.NumPartitions)
			}

			if kp.ElementsPerPart != tc.elementsPerPart {
				t.Errorf("Expected ElementsPerPart=%d, got %d",
					tc.elementsPerPart, kp.ElementsPerPart)
			}
		})
	}
}

// Test 2: Static Matrix Code Generation (Build Systematically)
func TestKernelProgram_TypeDefaults(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	// Test 1: Both unset (should get defaults)
	kp1 := NewKernelProgram(device, Config{
		Order:           1,
		ElementsPerPart: 1,
	})
	defer kp1.Free()

	if kp1.FloatType != Float64 {
		t.Errorf("Test 1: Expected FloatType=Float64(%d), got %d", Float64, kp1.FloatType)
	}
	if kp1.IntType != Int64 {
		t.Errorf("Test 1: Expected IntType=Int64(%d), got %d", Int64, kp1.IntType)
	}

	// Test 2: Explicit Float32 and Int32
	kp2 := NewKernelProgram(device, Config{
		Order:           1,
		ElementsPerPart: 1,
		FloatType:       Float32,
		IntType:         Int32,
	})
	defer kp2.Free()

	if kp2.FloatType != Float32 {
		t.Errorf("Test 2: Expected FloatType=Float32(%d), got %d", Float32, kp2.FloatType)
	}
	if kp2.IntType != Int32 {
		t.Errorf("Test 2: Expected IntType=Int32(%d), got %d", Int32, kp2.IntType)
	}

	// Test 3: Explicit Float64 and Int64
	kp3 := NewKernelProgram(device, Config{
		Order:           1,
		ElementsPerPart: 1,
		FloatType:       Float64,
		IntType:         Int64,
	})
	defer kp3.Free()

	if kp3.FloatType != Float64 {
		t.Errorf("Test 3: Expected FloatType=Float64(%d), got %d", Float64, kp3.FloatType)
	}
	if kp3.IntType != Int64 {
		t.Errorf("Test 3: Expected IntType=Int64(%d), got %d", Int64, kp3.IntType)
	}
}

func TestKernelProgram_StaticMatrixGeneration(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	kp := NewKernelProgram(device, Config{
		Order:           2,
		ElementsPerPart: 5,
		FloatType:       Float64,
	})
	defer kp.Free()

	// Create a simple 3x3 identity matrix
	identity := NewMatrix(3, 3, []float64{
		1.0, 0.0, 0.0,
		0.0, 1.0, 0.0,
		0.0, 0.0, 1.0,
	})

	kp.AddStaticMatrix("TestMatrix", identity)
	preamble := kp.GenerateKernelMain()

	// Check basic structure
	if !strings.Contains(preamble, "const double TestMatrix[3][3]") {
		t.Error("Preamble missing TestMatrix declaration")
	}

	// Check type definitions
	if !strings.Contains(preamble, "typedef double real_t") {
		t.Error("Preamble missing real_t typedef")
	}

	// Check constants
	if !strings.Contains(preamble, "#define NP") {
		t.Error("Preamble missing NP definition")
	}

	// Check partition-aware constants
	if !strings.Contains(preamble, "#define NPART") {
		t.Error("Preamble missing NPART definition")
	}
	if !strings.Contains(preamble, "#define K ") {
		t.Error("Preamble missing K (ElementsPerPart) definition")
	}

	// Check indexing macros
	if !strings.Contains(preamble, "#define NODE_IDX") {
		t.Error("Preamble missing NODE_IDX macro")
	}
}

// Test partition-aware memory allocation
func TestKernelProgram_PartitionMemoryAllocation(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	testCases := []struct {
		name            string
		numPartitions   int
		elementsPerPart int
		order           int
	}{
		{"SinglePartition", 1, 100, 3},
		{"TwoPartitions", 2, 50, 3},
		{"EightPartitions", 8, 25, 3},
		{"ManyPartitions", 64, 10, 3},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			kp := NewKernelProgram(device, Config{
				Order:           tc.order,
				NumPartitions:   tc.numPartitions,
				ElementsPerPart: tc.elementsPerPart,
				FloatType:       Float64,
			})
			defer kp.Free()

			err := kp.AllocateKernelMemory()
			if err != nil {
				t.Fatalf("Failed to allocate memory: %v", err)
			}

			// Verify total allocated size
			expectedNodes := tc.numPartitions * tc.elementsPerPart * kp.Np

			// Test that memory is correctly sized by copying data
			testData := make([]float64, expectedNodes)
			for i := range testData {
				testData[i] = float64(i)
			}

			kp.GetMemory("U").CopyFrom(unsafe.Pointer(&testData[0]),
				int64(expectedNodes*8))

			// Copy back and verify
			result := make([]float64, expectedNodes)
			kp.GetMemory("U").CopyTo(unsafe.Pointer(&result[0]),
				int64(expectedNodes*8))

			for i := range result {
				if result[i] != testData[i] {
					t.Errorf("Memory test failed at index %d: expected %f, got %f",
						i, testData[i], result[i])
					break
				}
			}
		})
	}
}

// Test 3: Index-Based Kernel (Specific Property Testing)
func TestKernelProgram_SimpleIndexKernel(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	kp := NewKernelProgram(device, Config{
		Order:           1,
		NumPartitions:   1,
		ElementsPerPart: 10,
		FloatType:       Float64,
	})
	defer kp.Free()

	// Generate preamble
	kp.GenerateKernelMain()

	// Create simple kernel that uses indices to set pattern
	kernelSource := `
	@kernel void indexPattern(
		const int N,
		const int_t* indices,
		real_t* data
	) {
		// Single outer loop iteration
		for (int block = 0; block < 1; ++block; @outer(0)) {
			// Inner loop does the actual work
			for (int i = 0; i < N; ++i; @inner(0)) {
				int idx = indices[i];
				real_t value = i * 2.0;
				data[idx] = value;
			}
		}
	}
	`

	_, err := kp.BuildKernel(kernelSource, "indexPattern")
	if err != nil {
		t.Fatalf("Failed to build kernel: %v", err)
	}

	// Create test data
	N := 10
	indices := []int64{9, 8, 7, 6, 5, 4, 3, 2, 1, 0} // Reverse order
	data := make([]float64, N)

	// Allocate device memory
	indicesMem := device.Malloc(int64(N*8), unsafe.Pointer(&indices[0]), nil)
	dataMem := device.Malloc(int64(N*8), unsafe.Pointer(&data[0]), nil)
	defer indicesMem.Free()
	defer dataMem.Free()

	// Execute kernel (note: RunKernel will set partition dims, but this kernel ignores them)
	err = kp.RunKernel("indexPattern", N, indicesMem, dataMem)
	if err != nil {
		t.Fatalf("Kernel execution failed: %v", err)
	}

	// Copy back and validate
	dataMem.CopyTo(unsafe.Pointer(&data[0]), int64(N*8))

	// Expected pattern: data[9]=0, data[8]=2, data[7]=4, etc.
	for i := 0; i < N; i++ {
		expected := float64((N - 1 - i) * 2)
		if math.Abs(data[i]-expected) > 1e-10 {
			t.Errorf("data[%d]: expected %.1f, got %.1f",
				i, expected, data[i])
		}
	}
}

// Test partition indexing macros
func TestKernelProgram_PartitionIndexing(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	kp := NewKernelProgram(device, Config{
		Order:           2,
		NumPartitions:   4,
		ElementsPerPart: 10,
		FloatType:       Float64,
	})
	defer kp.Free()

	preamble := kp.GenerateKernelMain()

	// Check for partition-aware constants
	if !strings.Contains(preamble, "#define NPART 4") {
		t.Error("Preamble missing NPART definition with correct value")
	}
	if !strings.Contains(preamble, "#define K 10") {
		t.Error("Preamble missing K (ElementsPerPart) definition with correct value")
	}

	// Test kernel using partition indexing
	kernelSource := `
	@kernel void partitionTest(
		const int Npart,
		const int K,
		real_t* data
	) {
		for (int part = 0; part < Npart; ++part; @outer(0)) {
			for (int elem = 0; elem < K; ++elem) {
				for (int node = 0; node < NP; ++node; @inner(0)) {
					int idx = NODE_IDX(part, elem, node);
					data[idx] = part * 1000.0 + elem * 10.0 + node;
				}
			}
		}
	}
	`

	_, err := kp.BuildKernel(kernelSource, "partitionTest")
	if err != nil {
		t.Fatalf("Failed to build partition kernel: %v", err)
	}

	// Allocate and test
	err = kp.AllocateKernelMemory()
	if err != nil {
		t.Fatalf("Failed to allocate memory: %v", err)
	}

	// Execute with partition dimensions
	err = kp.RunKernel("partitionTest",
		kp.NumPartitions, kp.ElementsPerPart, kp.GetMemory("U"))
	if err != nil {
		t.Fatalf("Kernel execution failed: %v", err)
	}

	// Verify partition-based pattern
	totalNodes := kp.NumPartitions * kp.ElementsPerPart * kp.Np
	result := make([]float64, totalNodes)
	kp.GetMemory("U").CopyTo(unsafe.Pointer(&result[0]), int64(totalNodes*8))

	// Check pattern
	for part := 0; part < kp.NumPartitions; part++ {
		for elem := 0; elem < kp.ElementsPerPart; elem++ {
			for node := 0; node < kp.Np; node++ {
				idx := part*kp.ElementsPerPart*kp.Np + elem*kp.Np + node
				expected := float64(part*1000 + elem*10 + node)
				if math.Abs(result[idx]-expected) > 1e-10 {
					t.Errorf("Partition %d, elem %d, node %d: expected %.0f, got %.0f",
						part, elem, node, expected, result[idx])
				}
			}
		}
	}
}

// Test 4: Matrix-Vector Product (Mathematical Exactness)
func TestKernelProgram_MatrixVectorProduct(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	// Test progressive orders
	for _, order := range []int{1, 2, 3} {
		t.Run(fmt.Sprintf("Order_%d", order), func(t *testing.T) {
			np := (order + 1) * (order + 2) * (order + 3) / 6

			kp := NewKernelProgram(device, Config{
				Order:           order,
				NumPartitions:   1,
				ElementsPerPart: 1,
				FloatType:       Float64,
			})
			defer kp.Free()

			// Create test differentiation matrix
			Dr := createTestDrMatrix(order, np)
			kp.AddStaticMatrix("Dr", Dr)
			kp.GenerateKernelMain()

			// Partition-aware kernel that applies Dr
			kernelSource := `
			@kernel void applyDr(
				const int Npart,
				const int K,
				const real_t* U,
				real_t* dU
			) {
				// Partition parallel execution
				for (int part = 0; part < Npart; ++part; @outer(0)) {
					// Apply Dr within this partition
					for (int node = 0; node < NP; ++node; @inner(0)) {
						// For this simple test, just apply within single partition
						if (part == 0) {
							matMul_Dr_Large(U, dU, K);
						}
					}
				}
			}
			`

			_, err := kp.BuildKernel(kernelSource, "applyDr")
			if err != nil {
				t.Fatalf("Failed to build kernel: %v", err)
			}

			// Test with polynomial that Dr should differentiate exactly
			U := createPolynomialData(order, np)
			expected := computeExpectedDerivative(Dr, U, np)

			// Allocate and copy
			err = kp.AllocateKernelMemory()
			if err != nil {
				t.Fatalf("Failed to allocate memory: %v", err)
			}

			kp.GetMemory("U").CopyFrom(unsafe.Pointer(&U[0]), int64(len(U)*8))

			// Execute
			err = kp.RunKernel("applyDr",
				kp.NumPartitions,
				kp.ElementsPerPart,
				kp.GetMemory("U"),
				kp.GetMemory("RHS"))
			if err != nil {
				t.Fatalf("Kernel execution failed: %v", err)
			}

			// Validate
			result := make([]float64, np)
			kp.GetMemory("RHS").CopyTo(unsafe.Pointer(&result[0]), int64(len(result)*8))

			for i := 0; i < np; i++ {
				if math.Abs(result[i]-expected[i]) > 1e-10 {
					t.Errorf("Order %d, node %d: expected %.6e, got %.6e",
						order, i, expected[i], result[i])
				}
			}
		})
	}
}

// Test 5: Chained Operations with Known Pattern (Complex Validation)
func TestKernelProgram_ChainedOperations(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	kp := NewKernelProgram(device, Config{
		Order:           2,
		NumPartitions:   1,
		ElementsPerPart: 5,
		FloatType:       Float64,
	})
	defer kp.Free()

	// Add multiple matrices
	kp.AddStaticMatrix("A", createPatternMatrix(kp.Np, 2.0))
	kp.AddStaticMatrix("B", createPatternMatrix(kp.Np, 3.0))
	kp.GenerateKernelMain()

	// Complex kernel with index indirection and chained ops
	kernelSource := `
@kernel void chainedPattern(
    const int K,
    const int_t* permutation,
    const real_t* input,
    real_t* temp1,
    real_t* temp2,
    real_t* output
) {
    // Single partition execution for this test
    for (int part = 0; part < 1; ++part; @outer(0)) {
        // Inner work over nodes
        for (int node = 0; node < NP; ++node; @inner(0)) {
            // Step 1: Apply matrix A
            matMul_A_Large(input, temp1, K);
            
            // Step 2: Permute using indices
            for (int i = 0; i < NP*K; ++i) {
                temp2[permutation[i]] = temp1[i];
            }
            
            // Step 3: Apply matrix B
            matMul_B_Large(temp2, output, K);
            
            // Step 4: Add signature pattern
            for (int i = 0; i < NP*K; ++i) {
                int mod_val = i % NP;
                output[i] += (real_t)mod_val;
            }
        }
    }
}
`

	_, err := kp.BuildKernel(kernelSource, "chainedPattern")
	if err != nil {
		t.Fatalf("Failed to build kernel: %v", err)
	}

	// Create permutation that reverses within each element
	perm := createElementWiseReversePermutation(kp.Np, kp.ElementsPerPart)

	// Create input with known pattern
	input := createCheckerboardPattern(kp.Np, kp.ElementsPerPart)

	// Allocate all memory
	err = kp.AllocateKernelMemory()
	if err != nil {
		t.Fatalf("Failed to allocate memory: %v", err)
	}

	// Additional allocations
	size := kp.Np * kp.ElementsPerPart
	permMem := device.Malloc(int64(size*8), unsafe.Pointer(&perm[0]), nil)
	temp1Mem := device.Malloc(int64(size*8), nil, nil)
	temp2Mem := device.Malloc(int64(size*8), nil, nil)
	defer permMem.Free()
	defer temp1Mem.Free()
	defer temp2Mem.Free()

	kp.GetMemory("U").CopyFrom(unsafe.Pointer(&input[0]), int64(len(input)*8))

	// Execute
	err = kp.RunKernel("chainedPattern",
		kp.ElementsPerPart,
		permMem,
		kp.GetMemory("U"),
		temp1Mem,
		temp2Mem,
		kp.GetMemory("RHS"),
	)
	if err != nil {
		t.Fatalf("Kernel execution failed: %v", err)
	}

	// Validate complex pattern
	result := make([]float64, len(input))
	kp.GetMemory("RHS").CopyTo(unsafe.Pointer(&result[0]), int64(len(result)*8))

	// Compute expected result on host
	expected := computeChainedResult(input, perm,
		kp.StaticMatrices["A"], kp.StaticMatrices["B"],
		kp.Np, kp.ElementsPerPart)

	maxError := 0.0
	for i := range result {
		error := math.Abs(result[i] - expected[i])
		if error > maxError {
			maxError = error
		}
	}

	if maxError > 1e-10 {
		t.Errorf("Maximum error %.6e exceeds tolerance", maxError)
		// Print first few differences for debugging
		for i := 0; i < min(10, len(result)); i++ {
			t.Logf("result[%d] = %.6e, expected = %.6e, diff = %.6e",
				i, result[i], expected[i], result[i]-expected[i])
		}
	}
}

// Test 6: Memory Allocation Sizes (Detailed Coverage)
func TestKernelProgram_MemoryAllocation(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	// Test both Float32 and Float64
	configs := []struct {
		name      string
		floatType DataType
		intType   DataType
		floatSize int
		intSize   int
	}{
		{"Float64_Int64", Float64, Int64, 8, 8},
		{"Float32_Int32", Float32, Int32, 4, 4},
	}

	for _, cfg := range configs {
		t.Run(cfg.name, func(t *testing.T) {
			kp := NewKernelProgram(device, Config{
				Order:           3,
				NumPartitions:   2,
				ElementsPerPart: 50,
				FloatType:       cfg.floatType,
				IntType:         cfg.intType,
			})
			defer kp.Free()

			// Verify size methods
			if kp.GetFloatSize() != cfg.floatSize {
				t.Errorf("Expected float size %d, got %d",
					cfg.floatSize, kp.GetFloatSize())
			}
			if kp.GetIntSize() != cfg.intSize {
				t.Errorf("Expected int size %d, got %d",
					cfg.intSize, kp.GetIntSize())
			}

			// Allocate memory
			err := kp.AllocateKernelMemory()
			if err != nil {
				t.Fatalf("Failed to allocate memory: %v", err)
			}

			// Verify all allocations exist
			expectedArrays := []string{
				"U", "RHS",
				"rx", "ry", "rz", "sx", "sy", "sz", "tx", "ty", "tz",
				"faceM", "faceP", "faceTypes",
				"nx", "ny", "nz",
				"Fscale", "bcData",
			}

			for _, name := range expectedArrays {
				if kp.GetMemory(name) == nil {
					t.Errorf("Memory allocation '%s' not found", name)
				}
			}
		})
	}
}

// Test multi-partition kernel execution
func TestKernelProgram_MultiPartitionExecution(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	kp := NewKernelProgram(device, Config{
		Order:           1,
		NumPartitions:   8,
		ElementsPerPart: 4,
		FloatType:       Float64,
	})
	defer kp.Free()

	kp.GenerateKernelMain()

	// Kernel that demonstrates partition-parallel execution
	kernelSource := `
	@kernel void partitionSum(
		const int Npart,
		const int K,
		real_t* data
	) {
		for (int part = 0; part < Npart; ++part; @outer(0)) {
			for (int node = 0; node < NP; ++node; @inner(0)) {
				// Sum pattern: each partition adds its ID to all its elements
				for (int elem = 0; elem < K; ++elem) {
					int idx = NODE_IDX(part, elem, node);
					data[idx] = (real_t)part + (real_t)elem * 0.1 + (real_t)node * 0.01;
				}
			}
		}
	}
	`

	_, err := kp.BuildKernel(kernelSource, "partitionSum")
	if err != nil {
		t.Fatalf("Failed to build kernel: %v", err)
	}

	err = kp.AllocateKernelMemory()
	if err != nil {
		t.Fatalf("Failed to allocate memory: %v", err)
	}

	// Execute
	err = kp.RunKernel("partitionSum",
		kp.NumPartitions, kp.ElementsPerPart, kp.GetMemory("U"))
	if err != nil {
		t.Fatalf("Kernel execution failed: %v", err)
	}

	// Verify results
	totalNodes := kp.NumPartitions * kp.ElementsPerPart * kp.Np
	result := make([]float64, totalNodes)
	kp.GetMemory("U").CopyTo(unsafe.Pointer(&result[0]), int64(totalNodes*8))

	// Check each partition has correct values
	for part := 0; part < kp.NumPartitions; part++ {
		for elem := 0; elem < kp.ElementsPerPart; elem++ {
			for node := 0; node < kp.Np; node++ {
				idx := part*kp.ElementsPerPart*kp.Np + elem*kp.Np + node
				expected := float64(part) + float64(elem)*0.1 + float64(node)*0.01
				if math.Abs(result[idx]-expected) > 1e-10 {
					t.Errorf("Part %d, elem %d, node %d: expected %f, got %f",
						part, elem, node, expected, result[idx])
				}
			}
		}
	}
}

// Helper functions for test data generation
func createTestDevice(t *testing.T) *gocca.OCCADevice {
	device, err := gocca.NewDevice(`{"mode": "Serial"}`)
	if err != nil {
		t.Fatalf("Failed to create device: %v", err)
	}
	return device
}

func createTestDrMatrix(order, np int) Matrix {
	// Create a simple differentiation matrix for testing
	// This is not a real Dr matrix, just for testing matrix operations
	data := make([]float64, np*np)
	for i := 0; i < np; i++ {
		for j := 0; j < np; j++ {
			if i == j {
				data[i*np+j] = -1.0
			} else if j == i+1 {
				data[i*np+j] = 1.0
			}
		}
	}
	return NewMatrix(np, np, data)
}

func createPolynomialData(order, np int) []float64 {
	// Create polynomial test data
	data := make([]float64, np)
	for i := 0; i < np; i++ {
		x := float64(i) / float64(np-1)
		data[i] = x * x // Simple quadratic
	}
	return data
}

func computeExpectedDerivative(Dr Matrix, U []float64, np int) []float64 {
	// Compute Dr * U
	result := make([]float64, np)
	for i := 0; i < np; i++ {
		sum := 0.0
		for j := 0; j < np; j++ {
			sum += Dr.At(i, j) * U[j]
		}
		result[i] = sum
	}
	return result
}

func createPatternMatrix(np int, scale float64) Matrix {
	// Create a pattern matrix for testing
	data := make([]float64, np*np)
	for i := 0; i < np; i++ {
		for j := 0; j < np; j++ {
			data[i*np+j] = scale * float64(i+j) / float64(np)
		}
	}
	return NewMatrix(np, np, data)
}

func createElementWiseReversePermutation(np, numElements int) []int64 {
	perm := make([]int64, np*numElements)
	for elem := 0; elem < numElements; elem++ {
		for i := 0; i < np; i++ {
			srcIdx := elem*np + i
			dstIdx := elem*np + (np - 1 - i)
			perm[srcIdx] = int64(dstIdx)
		}
	}
	return perm
}

func createCheckerboardPattern(np, numElements int) []float64 {
	data := make([]float64, np*numElements)
	for i := range data {
		if i%2 == 0 {
			data[i] = 1.0
		} else {
			data[i] = -1.0
		}
	}
	return data
}

func computeChainedResult(input []float64, perm []int64, A, B Matrix, np, numElements int) []float64 {
	// Step 1: Apply matrix A
	temp1 := make([]float64, len(input))
	for elem := 0; elem < numElements; elem++ {
		for i := 0; i < np; i++ {
			sum := 0.0
			for j := 0; j < np; j++ {
				sum += A.At(i, j) * input[elem*np+j]
			}
			temp1[elem*np+i] = sum
		}
	}

	// Step 2: Permute
	temp2 := make([]float64, len(input))
	for i := range temp1 {
		temp2[perm[i]] = temp1[i]
	}

	// Step 3: Apply matrix B
	result := make([]float64, len(input))
	for elem := 0; elem < numElements; elem++ {
		for i := 0; i < np; i++ {
			sum := 0.0
			for j := 0; j < np; j++ {
				sum += B.At(i, j) * temp2[elem*np+j]
			}
			result[elem*np+i] = sum
		}
	}

	// Step 4: Add signature pattern
	for i := range result {
		result[i] += float64(i % np)
	}

	return result
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

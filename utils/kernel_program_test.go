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
				Order:       order,
				NumElements: 10,
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
		})
	}
}

// Test 2: Static Matrix Code Generation (Build Systematically)
func TestKernelProgram_TypeDefaults(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	// Test 1: Both unset (should get defaults)
	kp1 := NewKernelProgram(device, Config{
		Order:       1,
		NumElements: 1,
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
		Order:       1,
		NumElements: 1,
		FloatType:   Float32,
		IntType:     Int32,
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
		Order:       1,
		NumElements: 1,
		FloatType:   Float64,
		IntType:     Int64,
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
		Order:       2,
		NumElements: 5,
		FloatType:   Float64,
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

	// Verify preamble contains expected content
	expectedPatterns := []string{
		"const double TestMatrix[3][3]",
		"1.000000000000000e+00",
		"matMul_TestMatrix_Large",
		"typedef double real_t",
		"#define NP 10", // Order 2 has Np=10
	}

	for _, pattern := range expectedPatterns {
		if !strings.Contains(preamble, pattern) {
			t.Errorf("Preamble missing expected pattern: %s", pattern)
		}
	}

	// Test Float32 generation
	kp32 := NewKernelProgram(device, Config{
		Order:       1,
		NumElements: 5,
		FloatType:   Float32,
		IntType:     Int32,
	})
	defer kp32.Free()

	// Verify types were set correctly
	if kp32.FloatType != Float32 {
		t.Errorf("Expected FloatType=Float32, got %v", kp32.FloatType)
	}
	if kp32.IntType != Int32 {
		t.Errorf("Expected IntType=Int32, got %v", kp32.IntType)
	}

	small := NewMatrix(2, 2, []float64{1.5, 2.5, 3.5, 4.5})
	kp32.AddStaticMatrix("SmallMat", small)
	preamble32 := kp32.GenerateKernelMain()

	// Debug: print preamble if test fails
	if len(preamble32) == 0 {
		t.Error("Float32 preamble is empty")
	}

	float32Patterns := []string{
		"__constant__ float SmallMat[2][2]",
		"typedef float real_t",
		"typedef int int_t",
		"1.5000000e+00f", // %.7ef format gives 7 digits after decimal
	}

	for _, pattern := range float32Patterns {
		if !strings.Contains(preamble32, pattern) {
			t.Errorf("Float32 preamble missing pattern: %s", pattern)
		}
	}
}

// Test 3: Simple Index-Based Kernel (Specific Property Testing)
func TestKernelProgram_SimpleIndexKernel(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	kp := NewKernelProgram(device, Config{
		Order:       1,
		NumElements: 10,
		FloatType:   Float64,
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

	// Execute kernel
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

// Test 4: Matrix-Vector Product (Mathematical Exactness)
func TestKernelProgram_MatrixVectorProduct(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	// Test progressive orders
	for _, order := range []int{1, 2, 3} {
		t.Run(fmt.Sprintf("Order_%d", order), func(t *testing.T) {
			np := (order + 1) * (order + 2) * (order + 3) / 6

			kp := NewKernelProgram(device, Config{
				Order:       order,
				NumElements: 1,
				FloatType:   Float64,
			})
			defer kp.Free()

			// Create test differentiation matrix
			Dr := createTestDrMatrix(order, np)
			kp.AddStaticMatrix("Dr", Dr)
			kp.GenerateKernelMain()

			// Kernel that applies Dr to data
			kernelSource := `
			@kernel void applyDr(
				const int K,
				const real_t* U,
				real_t* dU
			) {
				matMul_Dr_Large(U, dU, K);
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
			err = kp.RunKernel("applyDr", 1,
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
		Order:       2,
		NumElements: 5,
		FloatType:   Float64,
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
			output[i] += (real_t)(i % NP);
		}
	}
	`

	_, err := kp.BuildKernel(kernelSource, "chainedPattern")
	if err != nil {
		t.Fatalf("Failed to build kernel: %v", err)
	}

	// Create permutation that reverses within each element
	perm := createElementWiseReversePermutation(kp.Np, kp.NumElements)

	// Create input with known pattern
	input := createCheckerboardPattern(kp.Np, kp.NumElements)

	// Allocate all memory
	err = kp.AllocateKernelMemory()
	if err != nil {
		t.Fatalf("Failed to allocate memory: %v", err)
	}

	// Additional allocations
	size := kp.Np * kp.NumElements
	permMem := device.Malloc(int64(size*8), unsafe.Pointer(&perm[0]), nil)
	temp1Mem := device.Malloc(int64(size*8), nil, nil)
	temp2Mem := device.Malloc(int64(size*8), nil, nil)
	defer permMem.Free()
	defer temp1Mem.Free()
	defer temp2Mem.Free()

	kp.GetMemory("U").CopyFrom(unsafe.Pointer(&input[0]), int64(len(input)*8))

	// Execute
	err = kp.RunKernel("chainedPattern",
		kp.NumElements,
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
		kp.Np, kp.NumElements)

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
				Order:       2,
				NumElements: 10,
				FloatType:   cfg.floatType,
				IntType:     cfg.intType,
			})
			defer kp.Free()

			err := kp.AllocateKernelMemory()
			if err != nil {
				t.Fatalf("Failed to allocate memory: %v", err)
			}

			// Verify size methods
			if kp.GetFloatSize() != cfg.floatSize {
				t.Errorf("Expected float size %d, got %d",
					cfg.floatSize, kp.GetFloatSize())
			}
			if kp.GetIntSize() != cfg.intSize {
				t.Errorf("Expected int size %d, got %d",
					cfg.intSize, kp.GetIntSize())
			}

			// Verify key allocations exist
			requiredMemory := []string{
				"U", "RHS", "rx", "ry", "rz", "sx", "sy", "sz", "tx", "ty", "tz",
				"faceM", "faceP", "faceTypes", "nx", "ny", "nz", "Fscale", "bcData",
			}

			for _, name := range requiredMemory {
				if kp.GetMemory(name) == nil {
					t.Errorf("Required memory '%s' not allocated", name)
				}
			}
		})
	}
}

// Helper Functions

func createTestDevice(t *testing.T) *gocca.OCCADevice {
	device, err := gocca.NewDevice(`{"mode": "Serial"}`)
	if err != nil {
		t.Fatalf("Failed to create device: %v", err)
	}
	return device
}

func createTestDrMatrix(order, np int) Matrix {
	// Create a simple differentiation-like matrix
	data := make([]float64, np*np)
	for i := 0; i < np; i++ {
		for j := 0; j < np; j++ {
			if i == j {
				data[i*np+j] = -1.0
			} else if j == (i+1)%np {
				data[i*np+j] = 1.0
			}
		}
	}
	return NewMatrix(np, np, data)
}

func createPatternMatrix(size int, scale float64) Matrix {
	// Create a matrix with a recognizable pattern
	data := make([]float64, size*size)
	for i := 0; i < size; i++ {
		for j := 0; j < size; j++ {
			// Diagonal dominant with pattern
			if i == j {
				data[i*size+j] = scale
			} else {
				data[i*size+j] = scale * 0.1 * float64(i-j) / float64(size)
			}
		}
	}
	return NewMatrix(size, size, data)
}

func createPolynomialData(order, np int) []float64 {
	// Create data representing a polynomial of degree <= order
	data := make([]float64, np)
	for i := 0; i < np; i++ {
		// Simple polynomial that's easy to verify
		x := float64(i) / float64(np-1)
		data[i] = x*x + 2*x + 1
	}
	return data
}

func computeExpectedDerivative(Dr Matrix, U []float64, np int) []float64 {
	// Compute Dr * U on host for validation
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

func createElementWiseReversePermutation(np, numElements int) []int64 {
	// Create permutation that reverses order within each element
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
	// Create a pattern that alternates values
	data := make([]float64, np*numElements)
	for i := range data {
		elem := i / np
		node := i % np
		if (elem+node)%2 == 0 {
			data[i] = 1.0
		} else {
			data[i] = -1.0
		}
	}
	return data
}

func computeChainedResult(input []float64, perm []int64,
	A, B Matrix, np, numElements int) []float64 {

	size := np * numElements

	// Step 1: Apply matrix A (blocked by element)
	temp1 := make([]float64, size)
	for elem := 0; elem < numElements; elem++ {
		offset := elem * np
		for i := 0; i < np; i++ {
			sum := 0.0
			for j := 0; j < np; j++ {
				sum += A.At(i, j) * input[offset+j]
			}
			temp1[offset+i] = sum
		}
	}

	// Step 2: Permute
	temp2 := make([]float64, size)
	for i := 0; i < size; i++ {
		temp2[perm[i]] = temp1[i]
	}

	// Step 3: Apply matrix B
	result := make([]float64, size)
	for elem := 0; elem < numElements; elem++ {
		offset := elem * np
		for i := 0; i < np; i++ {
			sum := 0.0
			for j := 0; j < np; j++ {
				sum += B.At(i, j) * temp2[offset+j]
			}
			result[offset+i] = sum
		}
	}

	// Step 4: Add signature pattern
	for i := 0; i < size; i++ {
		result[i] += float64(i % np)
	}

	return result
}

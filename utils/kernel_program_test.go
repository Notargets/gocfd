package utils

import (
	"fmt"
	"math"
	"strings"
	"testing"
	"unsafe"

	"github.com/notargets/gocca"
)

// Test 1: Creation with Variable Partition Sizes
func TestKernelProgram_VariablePartitions(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	testCases := []struct {
		name string
		k    []int
	}{
		{"SinglePartition", []int{100}},
		{"EqualPartitions", []int{50, 50}},
		{"VariablePartitions", []int{30, 70, 45, 55}},
		{"ManyPartitions", []int{10, 20, 15, 25, 30, 35, 40, 45}},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			kp := NewKernelProgram(device, Config{
				K:         tc.k,
				FloatType: Float64,
				IntType:   Int64,
			})
			defer kp.Free()

			if kp.NumPartitions != len(tc.k) {
				t.Errorf("Expected NumPartitions=%d, got %d", len(tc.k), kp.NumPartitions)
			}

			for i, kval := range tc.k {
				if kp.K[i] != kval {
					t.Errorf("K[%d]: expected %d, got %d", i, kval, kp.K[i])
				}
			}
		})
	}
}

// Test 2: Array Allocation with Alignment
func TestKernelProgram_ArrayAllocation(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	kp := NewKernelProgram(device, Config{
		K:         []int{10, 20, 15},
		FloatType: Float64,
		IntType:   Int64,
	})
	defer kp.Free()

	// Define arrays with different alignments
	specs := []ArraySpec{
		{Name: "U", Size: 45 * 8 * 20, Alignment: CacheLineAlign}, // 45 elements * 8 bytes * 20 nodes
		{Name: "RHS", Size: 45 * 8 * 20, Alignment: NoAlignment},
		{Name: "workspace", Size: 1024 * 8, Alignment: PageAlign},
	}

	err := kp.AllocateArrays(specs)
	if err != nil {
		t.Fatalf("Failed to allocate arrays: %v", err)
	}

	// Verify allocations
	for _, spec := range specs {
		// Check _global allocation
		if kp.GetMemory(spec.Name) == nil {
			t.Errorf("Array %s_global not allocated", spec.Name)
		}

		// Check _offsets allocation
		if kp.pooledMemory[spec.Name+"_offsets"] == nil {
			t.Errorf("Array %s_offsets not allocated", spec.Name)
		}

		// Verify array is tracked
		found := false
		for _, name := range kp.allocatedArrays {
			if name == spec.Name {
				found = true
				break
			}
		}
		if !found {
			t.Errorf("Array %s not tracked in allocatedArrays", spec.Name)
		}
	}

	// Verify K array is allocated
	if kp.pooledMemory["K"] == nil {
		t.Error("K array not allocated")
	}
}

// Test 3: Partition Access Macro Generation
func TestKernelProgram_PartitionMacros(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	kp := NewKernelProgram(device, Config{
		K: []int{5, 10, 7},
	})
	defer kp.Free()

	// Add static matrix
	kp.AddStaticMatrix("Dr", createTestMatrix(4, 4))

	// Allocate arrays
	specs := []ArraySpec{
		{Name: "U", Size: 22 * 8 * 4, Alignment: NoAlignment},
		{Name: "RHS", Size: 22 * 8 * 4, Alignment: NoAlignment},
	}
	kp.AllocateArrays(specs)

	preamble := kp.GeneratePreamble()

	// Check partition access macros
	if !strings.Contains(preamble, "#define U_PART(part) (U_global + U_offsets[part])") {
		t.Error("Missing U_PART macro")
	}
	if !strings.Contains(preamble, "#define RHS_PART(part) (RHS_global + RHS_offsets[part])") {
		t.Error("Missing RHS_PART macro")
	}

	// Check vectorizable matrix macro
	if !strings.Contains(preamble, "#define MATMUL_Dr(IN, OUT, K_VAL, NP)") {
		t.Error("Missing MATMUL_Dr macro")
	}

	// Check NPART constant
	if !strings.Contains(preamble, "#define NPART 3") {
		t.Error("Missing or incorrect NPART definition")
	}
}

// Test 4: Variable Partition Kernel Execution
func TestKernelProgram_VariablePartitionExecution(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	// Variable sized partitions
	k := []int{5, 8, 3, 10}
	np := 4 // nodes per element

	kp := NewKernelProgram(device, Config{
		K:         k,
		FloatType: Float64,
		IntType:   Int64,
	})
	defer kp.Free()

	// Calculate total size
	totalElements := 0
	for _, kval := range k {
		totalElements += kval
	}
	totalNodes := totalElements * np

	// Allocate arrays
	specs := []ArraySpec{
		{Name: "data", Size: int64(totalNodes * 8), Alignment: CacheLineAlign},
	}
	err := kp.AllocateArrays(specs)
	if err != nil {
		t.Fatalf("Failed to allocate arrays: %v", err)
	}

	kp.GeneratePreamble()

	// Kernel that fills each partition with its ID
	kernelSource := `
@kernel void fillPartitions(
    const int_t* K,
    real_t* data_global,
    const int_t* data_offsets
) {
    for (int part = 0; part < NPART; ++part; @outer) {
        for (int i = 0; i < 32; ++i; @inner) {
            // Get partition data pointer
            real_t* data = data_PART(part);
            int k_part = K[part];
            
            // Fill all elements in this partition
            for (int elem = 0; elem < k_part; ++elem) {
                for (int node = 0; node < 4; ++node) {
                    int idx = elem * 4 + node;
                    if (i == 0) {  // Only first thread writes
                        data[idx] = (real_t)part;
                    }
                }
            }
        }
    }
}
`

	_, err = kp.BuildKernel(kernelSource, "fillPartitions")
	if err != nil {
		t.Fatalf("Failed to build kernel: %v", err)
	}

	// Execute kernel
	err = kp.RunKernel("fillPartitions", "data")
	if err != nil {
		t.Fatalf("Kernel execution failed: %v", err)
	}

	// Verify results
	result := make([]float64, totalNodes)
	kp.GetMemory("data").CopyTo(unsafe.Pointer(&result[0]), int64(totalNodes*8))

	idx := 0
	for part := 0; part < len(k); part++ {
		for elem := 0; elem < k[part]; elem++ {
			for node := 0; node < np; node++ {
				expected := float64(part)
				if math.Abs(result[idx]-expected) > 1e-10 {
					t.Errorf("Part %d, elem %d, node %d: expected %f, got %f",
						part, elem, node, expected, result[idx])
				}
				idx++
			}
		}
	}
}

// Test 5: Matrix Operations with Variable K
func TestKernelProgram_MatrixOperations(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	k := []int{2, 3} // Variable elements per partition
	np := 4

	kp := NewKernelProgram(device, Config{
		K:         k,
		FloatType: Float64,
		IntType:   Int64,
	})
	defer kp.Free()

	// Add differentiation matrix
	Dr := createIdentityMatrix(np)
	kp.AddStaticMatrix("Dr", Dr)

	totalElements := 5
	totalNodes := totalElements * np

	// Allocate arrays
	specs := []ArraySpec{
		{Name: "U", Size: int64(totalNodes * 8), Alignment: NoAlignment},
		{Name: "dU", Size: int64(totalNodes * 8), Alignment: NoAlignment},
	}
	err := kp.AllocateArrays(specs)
	if err != nil {
		t.Fatalf("Failed to allocate arrays: %v", err)
	}

	kp.GeneratePreamble()

	// Kernel using vectorizable macro
	kernelSource := `
@kernel void applyDr(
    const int_t* K,
    const real_t* U_global,
    const int_t* U_offsets,
    real_t* dU_global,
    const int_t* dU_offsets
) {
    for (int part = 0; part < NPART; ++part; @outer) {
        for (int tid = 0; tid < 32; ++tid; @inner) {
            if (tid == 0) {  // Single thread does the work
                const real_t* U = U_PART(part);
                real_t* dU = dU_PART(part);
                int k_part = K[part];
                
                // Apply differentiation matrix using macro
                MATMUL_Dr(U, dU, k_part, 4);
            }
        }
    }
}
`

	_, err = kp.BuildKernel(kernelSource, "applyDr")
	if err != nil {
		t.Fatalf("Failed to build kernel: %v", err)
	}

	// Create input data
	input := make([]float64, totalNodes)
	for i := range input {
		input[i] = float64(i)
	}
	kp.GetMemory("U").CopyFrom(unsafe.Pointer(&input[0]), int64(totalNodes*8))

	// Execute
	err = kp.RunKernel("applyDr", "U", "dU")
	if err != nil {
		t.Fatalf("Kernel execution failed: %v", err)
	}

	// Verify (identity matrix should copy input)
	result := make([]float64, totalNodes)
	kp.GetMemory("dU").CopyTo(unsafe.Pointer(&result[0]), int64(totalNodes*8))

	for i := range result {
		if math.Abs(result[i]-input[i]) > 1e-10 {
			t.Errorf("Index %d: expected %f, got %f", i, input[i], result[i])
		}
	}
}

// Test 6: Incremental Partition Testing
func TestKernelProgram_IncrementalPartitions(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	// Test partition counts from 1 to 6
	for numParts := 1; numParts <= 6; numParts++ {
		t.Run(fmt.Sprintf("%dPartitions", numParts), func(t *testing.T) {
			// Create K array with variable sizes
			k := make([]int, numParts)
			totalElems := 24
			for i := 0; i < numParts; i++ {
				k[i] = totalElems / numParts
				if i < totalElems%numParts {
					k[i]++
				}
			}

			kp := NewKernelProgram(device, Config{
				K:         k,
				FloatType: Float64,
				IntType:   Int64,
			})
			defer kp.Free()

			// Verify total elements preserved
			sum := 0
			for _, kval := range kp.K {
				sum += kval
			}
			if sum != totalElems {
				t.Errorf("Total elements %d != %d", sum, totalElems)
			}

			// Test simple allocation
			specs := []ArraySpec{
				{Name: "test", Size: int64(totalElems * 8), Alignment: NoAlignment},
			}
			err := kp.AllocateArrays(specs)
			if err != nil {
				t.Errorf("Failed to allocate for %d partitions: %v", numParts, err)
			}
		})
	}
}

// Test 7: Type Configuration
func TestKernelProgram_TypeConfiguration(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	configs := []struct {
		name      string
		floatType DataType
		intType   DataType
		floatStr  string
		intStr    string
	}{
		{"Float64_Int64", Float64, Int64, "double", "long"},
		{"Float32_Int32", Float32, Int32, "float", "int"},
		{"Float32_Int64", Float32, Int64, "float", "long"},
		{"Float64_Int32", Float64, Int32, "double", "int"},
	}

	for _, cfg := range configs {
		t.Run(cfg.name, func(t *testing.T) {
			kp := NewKernelProgram(device, Config{
				K:         []int{10},
				FloatType: cfg.floatType,
				IntType:   cfg.intType,
			})
			defer kp.Free()

			preamble := kp.GeneratePreamble()

			// Check type definitions
			if !strings.Contains(preamble, fmt.Sprintf("typedef %s real_t", cfg.floatStr)) {
				t.Errorf("Missing or incorrect real_t typedef for %s", cfg.floatStr)
			}
			if !strings.Contains(preamble, fmt.Sprintf("typedef %s int_t", cfg.intStr)) {
				t.Errorf("Missing or incorrect int_t typedef for %s", cfg.intStr)
			}

			// Check size methods
			expectedFloatSize := 8
			if cfg.floatType == Float32 {
				expectedFloatSize = 4
			}
			if kp.GetFloatSize() != expectedFloatSize {
				t.Errorf("Expected float size %d, got %d", expectedFloatSize, kp.GetFloatSize())
			}

			expectedIntSize := 8
			if cfg.intType == Int32 {
				expectedIntSize = 4
			}
			if kp.GetIntSize() != expectedIntSize {
				t.Errorf("Expected int size %d, got %d", expectedIntSize, kp.GetIntSize())
			}
		})
	}
}

// Test 8: Error Handling
func TestKernelProgram_ErrorHandling(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	t.Run("EmptyK", func(t *testing.T) {
		defer func() {
			if r := recover(); r == nil {
				t.Error("Expected panic for empty K array")
			}
		}()
		NewKernelProgram(device, Config{K: []int{}})
	})

	t.Run("NilDevice", func(t *testing.T) {
		defer func() {
			if r := recover(); r == nil {
				t.Error("Expected panic for nil device")
			}
		}()
		NewKernelProgram(nil, Config{K: []int{10}})
	})

	t.Run("MismatchedNumPartitions", func(t *testing.T) {
		defer func() {
			if r := recover(); r == nil {
				t.Error("Expected panic for mismatched NumPartitions")
			}
		}()
		NewKernelProgram(device, Config{
			NumPartitions: 3,
			K:             []int{10, 20}, // Only 2 elements
		})
	})
}

// Helper functions
func createTestDevice(t *testing.T) *gocca.OCCADevice {
	device, err := gocca.NewDevice(`{"mode": "Serial"}`)
	if err != nil {
		t.Fatalf("Failed to create device: %v", err)
	}
	return device
}

func createTestMatrix(rows, cols int) Matrix {
	data := make([]float64, rows*cols)
	for i := range data {
		data[i] = float64(i)
	}
	return NewMatrix(rows, cols, data)
}

func createIdentityMatrix(n int) Matrix {
	data := make([]float64, n*n)
	for i := 0; i < n; i++ {
		data[i*n+i] = 1.0
	}
	return NewMatrix(n, n, data)
}

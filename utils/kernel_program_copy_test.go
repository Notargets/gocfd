package utils

import (
	"fmt"
	"math"
	"testing"
	"unsafe"
)

// ============================================================================
// Section 1: Fundamental Tests - Basic Operations
// Following Unit Testing Principle: Start with fundamentals
// ============================================================================

// Test 1.1: GetArrayType basic functionality
func TestGetArrayType_BasicFunctionality(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	kp := NewKernelProgram(device, Config{
		K:         []int{10},
		FloatType: Float64,
		IntType:   Int32,
	})
	defer kp.Free()

	// Allocate arrays with different types
	specs := []ArraySpec{
		{Name: "float64_array", Size: 10 * 8, DataType: Float64, Alignment: NoAlignment},
		{Name: "float32_array", Size: 10 * 4, DataType: Float32, Alignment: NoAlignment},
		{Name: "int32_array", Size: 10 * 4, DataType: Int32, Alignment: NoAlignment},
		{Name: "int64_array", Size: 10 * 8, DataType: Int64, Alignment: NoAlignment},
	}

	err := kp.AllocateArrays(specs)
	if err != nil {
		t.Fatalf("Failed to allocate arrays: %v", err)
	}

	// Test each array type
	testCases := []struct {
		name     string
		expected DataType
	}{
		{"float64_array", Float64},
		{"float32_array", Float32},
		{"int32_array", Int32},
		{"int64_array", Int64},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			dataType, err := kp.GetArrayType(tc.name)
			if err != nil {
				t.Errorf("GetArrayType failed: %v", err)
			}
			if dataType != tc.expected {
				t.Errorf("Expected type %v, got %v", tc.expected, dataType)
			}
		})
	}
}

// Test 1.2: GetArrayType error cases
func TestGetArrayType_ErrorCases(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	kp := NewKernelProgram(device, Config{K: []int{10}})
	defer kp.Free()

	// Test non-existent array
	_, err := kp.GetArrayType("non_existent")
	if err == nil {
		t.Error("Expected error for non-existent array")
	}
}

// Test 1.3: GetArrayLogicalSize basic functionality
func TestGetArrayLogicalSize_SinglePartition(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	kp := NewKernelProgram(device, Config{
		K: []int{5},
	})
	defer kp.Free()

	// Test with different element counts per value
	testCases := []struct {
		name             string
		size             int64
		elementsPerValue int
		expectedSize     int
	}{
		{"scalar_array", 5 * 8, 1, 5},         // 5 elements * 1 value each
		{"vector3_array", 5 * 3 * 8, 3, 15},   // 5 elements * 3 values each
		{"matrix4_array", 5 * 16 * 8, 16, 80}, // 5 elements * 16 values each
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			spec := ArraySpec{
				Name:      tc.name,
				Size:      tc.size,
				DataType:  Float64,
				Alignment: NoAlignment,
			}

			// Manually set elementsPerValue in metadata after allocation
			err := kp.AllocateArrays([]ArraySpec{spec})
			if err != nil {
				t.Fatalf("Failed to allocate: %v", err)
			}

			// Update metadata to reflect elementsPerValue
			if metadata, exists := kp.arrayMetadata[tc.name]; exists {
				metadata.elementsPerValue = tc.elementsPerValue
				kp.arrayMetadata[tc.name] = metadata
			}

			size, err := kp.GetArrayLogicalSize(tc.name)
			if err != nil {
				t.Errorf("GetArrayLogicalSize failed: %v", err)
			}
			if size != tc.expectedSize {
				t.Errorf("Expected size %d, got %d", tc.expectedSize, size)
			}
		})
	}
}

// ============================================================================
// Section 2: Progressive Complexity - Multiple Partitions
// Following Unit Testing Principle: Build systematically
// ============================================================================

// Test 2.1: GetArrayLogicalSize with multiple partitions
func TestGetArrayLogicalSize_MultiplePartitions(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	// Test incrementally: 1, 2, 3, ... partitions
	testCases := []struct {
		name string
		k    []int
	}{
		{"2_partitions", []int{3, 4}},
		{"3_partitions", []int{2, 3, 4}},
		{"4_partitions", []int{1, 2, 3, 4}},
		{"5_partitions", []int{1, 1, 1, 1, 1}},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			kp := NewKernelProgram(device, Config{K: tc.k})
			defer kp.Free()

			totalK := 0
			for _, k := range tc.k {
				totalK += k
			}

			spec := ArraySpec{
				Name:      "test_array",
				Size:      int64(totalK * 8),
				DataType:  Float64,
				Alignment: NoAlignment,
			}

			err := kp.AllocateArrays([]ArraySpec{spec})
			if err != nil {
				t.Fatalf("Failed to allocate: %v", err)
			}

			size, err := kp.GetArrayLogicalSize("test_array")
			if err != nil {
				t.Errorf("GetArrayLogicalSize failed: %v", err)
			}

			// Mathematical property: size = sum(K[i]) * elementsPerValue
			if size != totalK {
				t.Errorf("Expected size %d, got %d", totalK, size)
			}
		})
	}
}

// Test 2.2: CopyArrayToHost single partition
func TestCopyArrayToHost_SinglePartition(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	kp := NewKernelProgram(device, Config{
		K:         []int{5},
		FloatType: Float64,
	})
	defer kp.Free()

	// Allocate and initialize array
	spec := ArraySpec{
		Name:      "data",
		Size:      5 * 8,
		DataType:  Float64,
		Alignment: NoAlignment,
	}
	err := kp.AllocateArrays([]ArraySpec{spec})
	if err != nil {
		t.Fatalf("Failed to allocate: %v", err)
	}

	// Write test data
	testData := []float64{1.0, 2.0, 3.0, 4.0, 5.0}
	mem := kp.GetMemory("data")
	mem.CopyFrom(unsafe.Pointer(&testData[0]), int64(5*8))

	// Copy back using generic method
	result, err := CopyArrayToHost[float64](kp, "data")
	if err != nil {
		t.Fatalf("CopyArrayToHost failed: %v", err)
	}

	// Verify
	if len(result) != 5 {
		t.Errorf("Expected 5 elements, got %d", len(result))
	}

	for i := 0; i < 5; i++ {
		if math.Abs(result[i]-testData[i]) > 1e-10 {
			t.Errorf("Element %d: expected %f, got %f", i, testData[i], result[i])
		}
	}
}

// ============================================================================
// Section 3: Type Safety and Verification
// Following Unit Testing Principle: Specific property testing
// ============================================================================

// Test 3.1: CopyArrayToHost type verification
func TestCopyArrayToHost_TypeVerification(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	kp := NewKernelProgram(device, Config{K: []int{10}})
	defer kp.Free()

	// Allocate Float64 array
	spec := ArraySpec{
		Name:      "float64_data",
		Size:      10 * 8,
		DataType:  Float64,
		Alignment: NoAlignment,
	}
	err := kp.AllocateArrays([]ArraySpec{spec})
	if err != nil {
		t.Fatalf("Failed to allocate: %v", err)
	}

	// Try to copy as wrong type (float32)
	_, err = CopyArrayToHost[float32](kp, "float64_data")
	if err == nil {
		t.Error("Expected type mismatch error")
	}

	// Try to copy as int32
	_, err = CopyArrayToHost[int32](kp, "float64_data")
	if err == nil {
		t.Error("Expected type mismatch error")
	}

	// Correct type should work
	_, err = CopyArrayToHost[float64](kp, "float64_data")
	if err != nil {
		t.Errorf("Correct type failed: %v", err)
	}
}

// Test 3.2: All supported types
func TestCopyArrayToHost_AllTypes(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	kp := NewKernelProgram(device, Config{K: []int{5}})
	defer kp.Free()

	// Test all data types
	testCases := []struct {
		name     string
		dataType DataType
		size     int64
		copyFunc func() (interface{}, error)
	}{
		{
			name:     "float32",
			dataType: Float32,
			size:     5 * 4,
			copyFunc: func() (interface{}, error) {
				return CopyArrayToHost[float32](kp, "float32")
			},
		},
		{
			name:     "float64",
			dataType: Float64,
			size:     5 * 8,
			copyFunc: func() (interface{}, error) {
				return CopyArrayToHost[float64](kp, "float64")
			},
		},
		{
			name:     "int32",
			dataType: Int32,
			size:     5 * 4,
			copyFunc: func() (interface{}, error) {
				return CopyArrayToHost[int32](kp, "int32")
			},
		},
		{
			name:     "int64",
			dataType: Int64,
			size:     5 * 8,
			copyFunc: func() (interface{}, error) {
				return CopyArrayToHost[int64](kp, "int64")
			},
		},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			spec := ArraySpec{
				Name:      tc.name,
				Size:      tc.size,
				DataType:  tc.dataType,
				Alignment: NoAlignment,
			}
			err := kp.AllocateArrays([]ArraySpec{spec})
			if err != nil {
				t.Fatalf("Failed to allocate: %v", err)
			}

			result, err := tc.copyFunc()
			if err != nil {
				t.Errorf("Copy failed: %v", err)
			}
			if result == nil {
				t.Error("Result is nil")
			}
		})
	}
}

// ============================================================================
// Section 4: Alignment and Padding Removal
// Following Unit Testing Principle: Test specific mathematical properties
// ============================================================================

// Test 4.1: CopyArrayToHost removes padding correctly
func TestCopyArrayToHost_PaddingRemoval(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	// Use odd-sized partitions to ensure padding
	k := []int{3, 5, 7}
	totalElements := 15

	kp := NewKernelProgram(device, Config{
		K:         k,
		FloatType: Float64,
	})
	defer kp.Free()

	// Allocate with 64-byte alignment
	spec := ArraySpec{
		Name:      "aligned_data",
		Size:      int64(totalElements * 8),
		DataType:  Float64,
		Alignment: CacheLineAlign,
	}
	err := kp.AllocateArrays([]ArraySpec{spec})
	if err != nil {
		t.Fatalf("Failed to allocate: %v", err)
	}

	// Get offsets to understand padding
	offsetsMem := kp.pooledMemory["aligned_data_offsets"]
	numOffsets := len(k) + 1
	offsets := make([]int64, numOffsets)
	offsetsMem.CopyTo(unsafe.Pointer(&offsets[0]), int64(numOffsets*8))

	// Write test pattern that includes padding areas
	_, totalSize := kp.calculateAlignedOffsetsAndSize(spec)
	paddedData := make([]float64, totalSize/8)

	// Fill with sentinel values
	for i := range paddedData {
		paddedData[i] = -999.0
	}

	// Fill actual data areas with sequential values
	dataValue := 0.0
	for part := 0; part < len(k); part++ {
		startIdx := offsets[part]
		for elem := 0; elem < k[part]; elem++ {
			paddedData[startIdx+int64(elem)] = dataValue
			dataValue++
		}
	}

	// Write to device
	mem := kp.GetMemory("aligned_data")
	mem.CopyFrom(unsafe.Pointer(&paddedData[0]), totalSize)

	// Copy back without padding
	result, err := CopyArrayToHost[float64](kp, "aligned_data")
	if err != nil {
		t.Fatalf("CopyArrayToHost failed: %v", err)
	}

	// Verify: should have exactly totalElements, no padding
	if len(result) != totalElements {
		t.Errorf("Expected %d elements, got %d", totalElements, len(result))
	}

	// Verify sequential values with no gaps
	for i := 0; i < totalElements; i++ {
		if math.Abs(result[i]-float64(i)) > 1e-10 {
			t.Errorf("Element %d: expected %f, got %f", i, float64(i), result[i])
		}
	}
}

// ============================================================================
// Section 5: CopyPartitionToHost Tests
// Following Unit Testing Principle: Incremental validation
// ============================================================================

// Test 5.1: CopyPartitionToHost single partition
func TestCopyPartitionToHost_SinglePartition(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	kp := NewKernelProgram(device, Config{K: []int{5}})
	defer kp.Free()

	spec := ArraySpec{
		Name:      "data",
		Size:      5 * 8,
		DataType:  Float64,
		Alignment: NoAlignment,
	}
	err := kp.AllocateArrays([]ArraySpec{spec})
	if err != nil {
		t.Fatalf("Failed to allocate: %v", err)
	}

	// Write test data
	testData := []float64{10, 20, 30, 40, 50}
	mem := kp.GetMemory("data")
	mem.CopyFrom(unsafe.Pointer(&testData[0]), int64(5*8))

	// Copy partition 0
	result, err := CopyPartitionToHost[float64](kp, "data", 0)
	if err != nil {
		t.Fatalf("CopyPartitionToHost failed: %v", err)
	}

	// Verify
	if len(result) != 5 {
		t.Errorf("Expected 5 elements, got %d", len(result))
	}

	for i := 0; i < 5; i++ {
		if math.Abs(result[i]-testData[i]) > 1e-10 {
			t.Errorf("Element %d: expected %f, got %f", i, testData[i], result[i])
		}
	}
}

// Test 5.2: CopyPartitionToHost multiple partitions incrementally
func TestCopyPartitionToHost_MultiplePartitions(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	k := []int{3, 4, 5}
	kp := NewKernelProgram(device, Config{
		K:         k,
		FloatType: Float64,
	})
	defer kp.Free()

	totalElements := 12
	spec := ArraySpec{
		Name:      "data",
		Size:      int64(totalElements * 8),
		DataType:  Float64,
		Alignment: NoAlignment,
	}
	err := kp.AllocateArrays([]ArraySpec{spec})
	if err != nil {
		t.Fatalf("Failed to allocate: %v", err)
	}

	// Write unique values to each partition
	fullData := make([]float64, totalElements)
	idx := 0
	for part, partK := range k {
		for elem := 0; elem < partK; elem++ {
			fullData[idx] = float64(part*100 + elem)
			idx++
		}
	}

	mem := kp.GetMemory("data")
	mem.CopyFrom(unsafe.Pointer(&fullData[0]), int64(totalElements*8))

	// Test each partition
	for partID, partK := range k {
		t.Run(fmt.Sprintf("partition_%d", partID), func(t *testing.T) {
			result, err := CopyPartitionToHost[float64](kp, "data", partID)
			if err != nil {
				t.Fatalf("CopyPartitionToHost failed: %v", err)
			}

			// Verify size
			if len(result) != partK {
				t.Errorf("Expected %d elements, got %d", partK, len(result))
			}

			// Verify values
			for elem := 0; elem < partK; elem++ {
				expected := float64(partID*100 + elem)
				if math.Abs(result[elem]-expected) > 1e-10 {
					t.Errorf("Partition %d, element %d: expected %f, got %f",
						partID, elem, expected, result[elem])
				}
			}
		})
	}
}

// Test 5.3: CopyPartitionToHost with alignment
func TestCopyPartitionToHost_WithAlignment(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	k := []int{3, 5, 7} // Odd sizes to test padding
	kp := NewKernelProgram(device, Config{
		K:         k,
		FloatType: Float64,
	})
	defer kp.Free()

	totalElements := 15
	spec := ArraySpec{
		Name:      "aligned_data",
		Size:      int64(totalElements * 8),
		DataType:  Float64,
		Alignment: CacheLineAlign,
	}
	err := kp.AllocateArrays([]ArraySpec{spec})
	if err != nil {
		t.Fatalf("Failed to allocate: %v", err)
	}

	// Get aligned memory layout info
	_, totalSize := kp.calculateAlignedOffsetsAndSize(spec)
	offsetsMem := kp.pooledMemory["aligned_data_offsets"]
	offsets := make([]int64, len(k)+1)
	offsetsMem.CopyTo(unsafe.Pointer(&offsets[0]), int64((len(k)+1)*8))

	// Create padded buffer with test data
	paddedBuffer := make([]float64, totalSize/8)

	// Fill with sentinel values
	for i := range paddedBuffer {
		paddedBuffer[i] = -1.0
	}

	// Write actual data
	for part, partK := range k {
		startIdx := offsets[part]
		for elem := 0; elem < partK; elem++ {
			paddedBuffer[startIdx+int64(elem)] = float64(part*1000 + elem)
		}
	}

	// Write to device
	mem := kp.GetMemory("aligned_data")
	mem.CopyFrom(unsafe.Pointer(&paddedBuffer[0]), totalSize)

	// Test each partition
	for partID, partK := range k {
		t.Run(fmt.Sprintf("aligned_partition_%d", partID), func(t *testing.T) {
			result, err := CopyPartitionToHost[float64](kp, "aligned_data", partID)
			if err != nil {
				t.Fatalf("CopyPartitionToHost failed: %v", err)
			}

			// Should get exactly partK elements, no padding
			if len(result) != partK {
				t.Errorf("Expected %d elements, got %d", partK, len(result))
			}

			// Verify correct values
			for elem := 0; elem < partK; elem++ {
				expected := float64(partID*1000 + elem)
				if math.Abs(result[elem]-expected) > 1e-10 {
					t.Errorf("Element %d: expected %f, got %f", elem, expected, result[elem])
				}
			}
		})
	}
}

// ============================================================================
// Section 6: Error Cases and Edge Conditions
// Following Unit Testing Principle: Detailed coverage of error conditions
// ============================================================================

// Test 6.1: CopyPartitionToHost error cases
func TestCopyPartitionToHost_ErrorCases(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	kp := NewKernelProgram(device, Config{K: []int{3, 4}})
	defer kp.Free()

	spec := ArraySpec{
		Name:      "data",
		Size:      7 * 8,
		DataType:  Float64,
		Alignment: NoAlignment,
	}
	err := kp.AllocateArrays([]ArraySpec{spec})
	if err != nil {
		t.Fatalf("Failed to allocate: %v", err)
	}

	// Test invalid partition ID (negative)
	_, err = CopyPartitionToHost[float64](kp, "data", -1)
	if err == nil {
		t.Error("Expected error for negative partition ID")
	}

	// Test invalid partition ID (too large)
	_, err = CopyPartitionToHost[float64](kp, "data", 2)
	if err == nil {
		t.Error("Expected error for partition ID >= NumPartitions")
	}

	// Test non-existent array
	_, err = CopyPartitionToHost[float64](kp, "non_existent", 0)
	if err == nil {
		t.Error("Expected error for non-existent array")
	}

	// Test type mismatch
	_, err = CopyPartitionToHost[int32](kp, "data", 0)
	if err == nil {
		t.Error("Expected type mismatch error")
	}
}

// Test 6.2: Edge case - empty partitions
func TestCopyMethods_EmptyPartitions(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	// Include zero-sized partitions
	k := []int{0, 3, 0, 2}
	kp := NewKernelProgram(device, Config{K: k})
	defer kp.Free()

	totalElements := 5
	spec := ArraySpec{
		Name:      "data",
		Size:      int64(totalElements * 8),
		DataType:  Float64,
		Alignment: NoAlignment,
	}
	err := kp.AllocateArrays([]ArraySpec{spec})
	if err != nil {
		t.Fatalf("Failed to allocate: %v", err)
	}

	// Test CopyPartitionToHost on empty partition
	result, err := CopyPartitionToHost[float64](kp, "data", 0)
	if err != nil {
		t.Fatalf("CopyPartitionToHost failed on empty partition: %v", err)
	}
	if len(result) != 0 {
		t.Errorf("Expected 0 elements for empty partition, got %d", len(result))
	}

	// Test GetArrayLogicalSize
	size, err := kp.GetArrayLogicalSize("data")
	if err != nil {
		t.Fatalf("GetArrayLogicalSize failed: %v", err)
	}
	if size != totalElements {
		t.Errorf("Expected size %d, got %d", totalElements, size)
	}
}

// ============================================================================
// Section 7: Integration Tests with Different Offset Types
// Following Unit Testing Principle: Real-world scenarios
// ============================================================================

// Test 7.1: Test with Int32 offsets
func TestCopyMethods_Int32Offsets(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	kp := NewKernelProgram(device, Config{
		K:       []int{100, 200, 150},
		IntType: Int32, // Force 32-bit offsets
	})
	defer kp.Free()

	totalElements := 450
	spec := ArraySpec{
		Name:      "large_data",
		Size:      int64(totalElements * 8),
		DataType:  Float64,
		Alignment: NoAlignment,
	}
	err := kp.AllocateArrays([]ArraySpec{spec})
	if err != nil {
		t.Fatalf("Failed to allocate: %v", err)
	}

	// Initialize with sequential values
	testData := make([]float64, totalElements)
	for i := 0; i < totalElements; i++ {
		testData[i] = float64(i)
	}
	mem := kp.GetMemory("large_data")
	mem.CopyFrom(unsafe.Pointer(&testData[0]), int64(totalElements*8))

	// Test CopyArrayToHost
	result, err := CopyArrayToHost[float64](kp, "large_data")
	if err != nil {
		t.Fatalf("CopyArrayToHost failed: %v", err)
	}

	// Verify all data
	for i := 0; i < totalElements; i++ {
		if math.Abs(result[i]-float64(i)) > 1e-10 {
			t.Errorf("Element %d: expected %f, got %f", i, float64(i), result[i])
		}
	}

	// Test partition copy with int32 offsets
	part1Result, err := CopyPartitionToHost[float64](kp, "large_data", 1)
	if err != nil {
		t.Fatalf("CopyPartitionToHost failed: %v", err)
	}
	if len(part1Result) != 200 {
		t.Errorf("Expected 200 elements in partition 1, got %d", len(part1Result))
	}
}

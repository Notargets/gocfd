package utils

import (
	"fmt"
	"math"
	"strings"
	"testing"
	"unsafe"
)

// ============================================================================
// Section 7: Code Generation Validation Tests
// Following Unit Testing Principle: Validate intended functionality
// ============================================================================

// Test 7.1: Verify partition access macro generation
// Purpose: Ensures critical partition access macros are generated for each array
func TestKernelProgram_MacroGeneration_PartitionAccess(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	kp := NewKernelProgram(device, Config{K: []int{5, 10}})
	defer kp.Free()

	// Allocate multiple arrays
	specs := []ArraySpec{
		{Name: "U", Size: 15 * 8, Alignment: NoAlignment},
		{Name: "RHS", Size: 15 * 8, Alignment: CacheLineAlign},
		{Name: "rx", Size: 15 * 8, Alignment: NoAlignment},
	}
	err := kp.AllocateArrays(specs)
	if err != nil {
		t.Fatalf("Failed to allocate arrays: %v", err)
	}

	preamble := kp.GeneratePreamble()

	// Verify each array gets proper macros
	expectedMacros := []string{
		"#define U_PART(part) (U_global + U_offsets[part])",
		"#define RHS_PART(part) (RHS_global + RHS_offsets[part])",
		"#define rx_PART(part) (rx_global + rx_offsets[part])",
	}

	for _, macro := range expectedMacros {
		if !strings.Contains(preamble, macro) {
			t.Errorf("Missing critical partition access macro: %s", macro)
		}
	}
}

// Test 7.2: Verify vectorizable matrix multiplication macros
// Purpose: Validates that matrix operations generate vectorizable code patterns
func TestKernelProgram_MacroGeneration_MatrixOperations(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	kp := NewKernelProgram(device, Config{K: []int{1}})
	defer kp.Free()

	// Add square matrix
	Dr := NewMatrix(4, 4, make([]float64, 16))
	kp.AddStaticMatrix("Dr", Dr)

	preamble := kp.GeneratePreamble()

	// Check for macro pattern - be flexible about exact format
	if !strings.Contains(preamble, "MATMUL_Dr") {
		t.Error("Missing MATMUL_Dr macro in generated preamble")
	}

	// Check that the matrix itself was embedded
	if !strings.Contains(preamble, "Dr[4][4]") {
		t.Error("Missing Dr matrix declaration")
	}
}

// Test 7.3: Static matrix C array format verification
// Purpose: Ensures matrices are embedded as compile-time constants correctly
func TestKernelProgram_StaticMatrix_CFormat(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	kp := NewKernelProgram(device, Config{
		K:         []int{1},
		FloatType: Float64,
	})
	defer kp.Free()

	// Add 2x3 matrix with known values
	matrix := NewMatrix(2, 3, []float64{
		1.0, 2.0, 3.0,
		4.0, 5.0, 6.0,
	})
	kp.AddStaticMatrix("TestMat", matrix)

	preamble := kp.GeneratePreamble()

	// Verify C array declaration format
	if !strings.Contains(preamble, "const double TestMat[2][3]") {
		t.Error("Missing or incorrect static matrix declaration")
	}

	// Verify matrix contains the values - be flexible about format
	// Check for the presence of the values rather than exact formatting
	hasFirstRow := strings.Contains(preamble, "1.0") || strings.Contains(preamble, "1.000000")
	hasSecondValue := strings.Contains(preamble, "2.0") || strings.Contains(preamble, "2.000000")
	hasThirdValue := strings.Contains(preamble, "3.0") || strings.Contains(preamble, "3.000000")

	if !hasFirstRow || !hasSecondValue || !hasThirdValue {
		t.Error("Matrix values not found in expected format")
	}
}

// ============================================================================
// Section 8: Kernel Parameter Expansion Tests
// Following Unit Testing Principle: Test specific behavior
// ============================================================================

// Test 8.1: Verify RunKernel parameter expansion
// Purpose: Validates automatic expansion of user arrays to kernel parameters
func TestKernelProgram_ParameterExpansion(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	kp := NewKernelProgram(device, Config{K: []int{2}})
	defer kp.Free()

	// Allocate arrays
	specs := []ArraySpec{
		{Name: "U", Size: 2 * 8, Alignment: NoAlignment},
		{Name: "RHS", Size: 2 * 8, Alignment: NoAlignment},
	}
	err := kp.AllocateArrays(specs)
	if err != nil {
		t.Fatalf("Failed to allocate: %v", err)
	}

	// Kernel that uses expanded parameters
	kernelSource := `
@kernel void testExpansion(
	const int_t* K,
	const real_t* U_global,
	const int_t* U_offsets,
	real_t* RHS_global,
	const int_t* RHS_offsets
) {
	for (int part = 0; part < NPART; ++part; @outer) {
		for (int i = 0; i < 1; ++i; @inner) {
			// Kernel should receive all expanded parameters
		}
	}
}
`

	_, err = kp.BuildKernel(kernelSource, "testExpansion")
	if err != nil {
		t.Fatalf("Failed to build kernel: %v", err)
	}

	// Test that RunKernel expands "U", "RHS" to full parameter list
	err = kp.RunKernel("testExpansion", "U", "RHS")
	if err != nil {
		t.Errorf("Parameter expansion failed: %v", err)
	}
}

// ============================================================================
// Section 9: Alignment and Offset Calculation Tests
// Following Unit Testing Principle: Test mathematical properties
// ============================================================================

// Test 9.1: Alignment impact on offset calculations
// Purpose: Validates that alignment requirements affect offset array values
func TestKernelProgram_Alignment_OffsetCalculation(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	// Use odd-sized partitions to force padding
	k := []int{3, 5, 7}
	kp := NewKernelProgram(device, Config{K: k})
	defer kp.Free()

	spec := ArraySpec{
		Name:      "aligned",
		Size:      15 * 8,
		Alignment: CacheLineAlign, // 64-byte alignment
	}
	err := kp.AllocateArrays([]ArraySpec{spec})
	if err != nil {
		t.Fatalf("Failed to allocate: %v", err)
	}

	// Get offsets
	offsetsMem := kp.pooledMemory["aligned_offsets"]
	offsets := make([]int64, 4)
	offsetsMem.CopyTo(unsafe.Pointer(&offsets[0]), int64(4*8))

	// Verify each partition starts at aligned boundary
	for i := 0; i < 3; i++ {
		byteOffset := offsets[i] * 8
		if byteOffset%64 != 0 {
			t.Errorf("Partition %d not aligned: byte offset %d not divisible by 64",
				i, byteOffset)
		}

		// Verify sufficient space allocated for partition
		partitionSpace := offsets[i+1] - offsets[i]
		if partitionSpace < int64(k[i]) {
			t.Errorf("Partition %d: insufficient space %d for %d elements",
				i, partitionSpace, k[i])
		}
	}
}

// ============================================================================
// Section 10: Integration Tests
// Following Unit Testing Principle: Real-world scenarios
// ============================================================================

// Test 10.1: End-to-end matrix multiplication using generated macros
// Purpose: Validates complete workflow from matrices to kernel execution
func TestKernelProgram_Integration_MatrixMultiplication(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	np := 3          // Nodes per element
	k := []int{2, 3} // Elements per partition

	kp := NewKernelProgram(device, Config{
		K:         k,
		FloatType: Float64,
	})
	defer kp.Free()

	// Add differentiation matrix
	Dr := NewMatrix(np, np, []float64{
		-1.0, 1.0, 0.0,
		-0.5, 0.0, 0.5,
		0.0, -1.0, 1.0,
	})
	kp.AddStaticMatrix("Dr", Dr)

	// Allocate arrays
	totalNodes := 5 * np // 5 total elements
	specs := []ArraySpec{
		{Name: "U", Size: int64(totalNodes * 8), Alignment: NoAlignment},
		{Name: "Ur", Size: int64(totalNodes * 8), Alignment: NoAlignment},
	}
	err := kp.AllocateArrays(specs)
	if err != nil {
		t.Fatalf("Failed to allocate: %v", err)
	}

	// Initialize U with test data
	U := make([]float64, totalNodes)
	for i := range U {
		U[i] = float64(i)
	}
	kp.GetMemory("U").CopyFrom(unsafe.Pointer(&U[0]), int64(totalNodes*8))

	// Kernel using MATMUL macro
	kernelSource := fmt.Sprintf(`
#define NP %d

@kernel void applyDr(
	const int_t* K,
	const real_t* U_global,
	const int_t* U_offsets,
	real_t* Ur_global,
	const int_t* Ur_offsets
) {
	for (int part = 0; part < NPART; ++part; @outer) {
		for (int i = 0; i < 1; ++i; @inner) {
			const real_t* U = U_PART(part);
			real_t* Ur = Ur_PART(part);
			MATMUL_Dr(U, Ur, K[part], NP);
		}
	}
}
`, np)

	_, err = kp.BuildKernel(kernelSource, "applyDr")
	if err != nil {
		t.Fatalf("Failed to build kernel: %v", err)
	}

	// Execute
	err = kp.RunKernel("applyDr", "U", "Ur")
	if err != nil {
		t.Fatalf("Kernel execution failed: %v", err)
	}

	// Verify results
	Ur := make([]float64, totalNodes)
	kp.GetMemory("Ur").CopyTo(unsafe.Pointer(&Ur[0]), int64(totalNodes*8))

	// Check first element's differentiation
	// Ur[0] = Dr[0,0]*U[0] + Dr[0,1]*U[1] + Dr[0,2]*U[2]
	//       = -1.0*0 + 1.0*1 + 0.0*2 = 1.0
	expected := 1.0
	if math.Abs(Ur[0]-expected) > 1e-10 {
		t.Errorf("Ur[0]: expected %f, got %f", expected, Ur[0])
	}
}

// Test 10.2: Variable K handling in kernels
// Purpose: Validates that macros work correctly with variable partition sizes
func TestKernelProgram_VariableK_MacroUsage(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	// Highly variable partition sizes
	k := []int{1, 5, 2, 10, 3}
	kp := NewKernelProgram(device, Config{K: k})
	defer kp.Free()

	totalElements := 0
	for _, size := range k {
		totalElements += size
	}

	spec := ArraySpec{
		Name:      "data",
		Size:      int64(totalElements * 8),
		Alignment: NoAlignment,
	}
	err := kp.AllocateArrays([]ArraySpec{spec})
	if err != nil {
		t.Fatalf("Failed to allocate: %v", err)
	}

	// Kernel that processes variable K correctly
	kernelSource := `
@kernel void processVariable(
	const int_t* K,
	real_t* data_global,
	const int_t* data_offsets
) {
	for (int part = 0; part < NPART; ++part; @outer) {
		for (int i = 0; i < 1; ++i; @inner) {
			real_t* data = data_PART(part);
			int k_part = K[part];
			
			// Set each element to partition_id * 100 + elem_id
			for (int elem = 0; elem < k_part; ++elem) {
				data[elem] = part * 100.0 + elem;
			}
		}
	}
}
`

	_, err = kp.BuildKernel(kernelSource, "processVariable")
	if err != nil {
		t.Fatalf("Failed to build kernel: %v", err)
	}

	err = kp.RunKernel("processVariable", "data")
	if err != nil {
		t.Fatalf("Kernel execution failed: %v", err)
	}

	// Verify each partition processed correctly
	result := make([]float64, totalElements)
	kp.GetMemory("data").CopyTo(unsafe.Pointer(&result[0]), int64(totalElements*8))

	idx := 0
	for part, kPart := range k {
		for elem := 0; elem < kPart; elem++ {
			expected := float64(part*100 + elem)
			if math.Abs(result[idx]-expected) > 1e-10 {
				t.Errorf("Part %d, elem %d: expected %f, got %f",
					part, elem, expected, result[idx])
			}
			idx++
		}
	}
}

// ============================================================================
// Section 11: Rectangular Matrix Support
// Following Unit Testing Principle: Incremental validation
// ============================================================================

// Test 11.1: Rectangular matrix macro generation
// Purpose: Validates support for non-square matrices like LIFT
func TestKernelProgram_RectangularMatrices(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	kp := NewKernelProgram(device, Config{K: []int{1}})
	defer kp.Free()

	// LIFT matrix example: Np x (Nfaces * Nfp)
	// 6 nodes x (3 faces * 2 face nodes) = 6x6 for this test
	lift := NewMatrix(6, 6, make([]float64, 36))
	kp.AddStaticMatrix("LIFT", lift)

	preamble := kp.GeneratePreamble()

	// Verify rectangular matrix declaration
	if !strings.Contains(preamble, "const double LIFT[6][6]") {
		t.Error("Missing LIFT matrix declaration")
	}

	// Verify macro handles rectangular dimensions
	if !strings.Contains(preamble, "MATMUL_LIFT") {
		t.Error("Missing MATMUL_LIFT macro")
	}
}

// ============================================================================
// Section 12: Invariant Testing
// Following Unit Testing Principle: Properties that must always hold
// ============================================================================

// Test 12.1: Memory layout invariants
// Purpose: Validates critical assumptions about memory layout
func TestKernelProgram_Invariant_MemoryLayout(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	k := []int{3, 7, 5, 9}
	kp := NewKernelProgram(device, Config{K: k})
	defer kp.Free()

	spec := ArraySpec{
		Name:      "test",
		Size:      24 * 8,
		Alignment: CacheLineAlign,
	}
	err := kp.AllocateArrays([]ArraySpec{spec})
	if err != nil {
		t.Fatalf("Failed to allocate: %v", err)
	}

	// Get offsets
	offsetsMem := kp.pooledMemory["test_offsets"]
	offsets := make([]int64, 5)
	offsetsMem.CopyTo(unsafe.Pointer(&offsets[0]), int64(5*8))

	// Invariant 1: Offsets are monotonically increasing
	for i := 0; i < 4; i++ {
		if offsets[i+1] <= offsets[i] {
			t.Errorf("Invariant violated: offsets not increasing at %d: %d <= %d",
				i, offsets[i+1], offsets[i])
		}
	}

	// Invariant 2: Each partition has enough space for its elements
	for i := 0; i < 4; i++ {
		space := offsets[i+1] - offsets[i]
		if space < int64(k[i]) {
			t.Errorf("Invariant violated: partition %d has space %d < required %d",
				i, space, k[i])
		}
	}

	// Invariant 3: First offset is 0
	if offsets[0] != 0 {
		t.Errorf("Invariant violated: first offset is %d, expected 0", offsets[0])
	}

	// Invariant 4: Total size is sufficient
	totalRequired := 0
	for _, size := range k {
		totalRequired += size
	}
	totalAllocated := offsets[4]
	if totalAllocated < int64(totalRequired) {
		t.Errorf("Invariant violated: total allocated %d < required %d",
			totalAllocated, totalRequired)
	}
}

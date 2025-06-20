package utils

import (
	"math"
	"strings"
	"testing"
	"unsafe"
)

// Integration test to ensure kernels compile with negative values
func TestKernelProgram_CompileWithNegativeValues(t *testing.T) {
	device := createTestDevice(t)
	defer device.Free()

	kp := NewKernelProgram(device, Config{
		Order:       1,
		NumElements: 4,
		FloatType:   Float64,
	})
	defer kp.Free()

	// Add a matrix with negative values
	testMatrix := NewMatrix(2, 2, []float64{1.0, -1.0, -1.0, 1.0})
	kp.AddStaticMatrix("SignMat", testMatrix)

	// Generate preamble
	preamble := kp.GenerateKernelMain()

	// Log the generated matrix for debugging
	t.Log("Generated SignMat declaration:")
	lines := strings.Split(preamble, "\n")
	matrixStarted := false
	for _, line := range lines {
		if strings.Contains(line, "const double SignMat[2][2]") {
			matrixStarted = true
		}
		if matrixStarted {
			t.Log(line)
			if strings.Contains(line, "};") {
				break
			}
		}
	}

	// Simple kernel that uses FIXED indices (OCCA limitation)
	kernelSource := `
	@kernel void useStaticMatrix(
		const int N,
		const real_t* input,
		real_t* output
	) {
		for (int elem = 0; elem < 1; ++elem; @outer(0)) {
			for (int idx = 0; idx < 1; ++idx; @inner(0)) {
				// Use fixed indices to avoid OCCA issues
				output[0] = SignMat[0][0] * input[0] + SignMat[0][1] * input[1];
				output[1] = SignMat[1][0] * input[0] + SignMat[1][1] * input[1];
			}
		}
	}
	`

	// This should compile successfully with the negative value fix
	_, err := kp.BuildKernel(kernelSource, "useStaticMatrix")
	if err != nil {
		t.Fatalf("Failed to build kernel with negative values: %v", err)
	}
	// Don't defer kernel.Free() - kp.Free() will handle it

	// Test execution
	N := 2
	input := []float64{2.0, 3.0}
	output := make([]float64, N)

	// Allocate device memory
	inputMem := device.Malloc(int64(N*8), unsafe.Pointer(&input[0]), nil)
	outputMem := device.Malloc(int64(N*8), nil, nil)
	defer inputMem.Free()
	defer outputMem.Free()

	// Execute kernel
	err = kp.RunKernel("useStaticMatrix", N, inputMem, outputMem)
	if err != nil {
		t.Fatalf("Kernel execution failed: %v", err)
	}

	// Copy back and validate
	outputMem.CopyTo(unsafe.Pointer(&output[0]), int64(N*8))

	// Expected: output[0] = 1*2 + (-1)*3 = -1
	//           output[1] = (-1)*2 + 1*3 = 1
	expected := []float64{-1.0, 1.0}
	for i := 0; i < N; i++ {
		if math.Abs(output[i]-expected[i]) > 1e-10 {
			t.Errorf("Output[%d] = %f, expected %f", i, output[i], expected[i])
		}
	}
}

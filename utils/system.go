package utils

import (
	"fmt"
	"math"
	"runtime"
)

func GetMemUsage() string {
	var m runtime.MemStats
	runtime.ReadMemStats(&m)
	// For info on each, see: https://golang.org/pkg/runtime/#MemStats
	bToMb := func(b uint64) uint64 {
		return b / 1024 / 1024
	}
	return fmt.Sprintf("Alloc = %v MiB TotalAlloc = %v MiB Sys = %v MiB NumGC = %v",
		bToMb(m.Alloc), bToMb(m.TotalAlloc), bToMb(m.Sys), m.NumGC)
}

func IsNanPanic(A any) {
	if IsNan(A) {
		panic("NAN found")
	}
}

func IsNan(A any) bool {
	switch v := A.(type) {
	case float64:
		return math.IsNaN(float64(v))
	case float32:
		return math.IsNaN(float64(v))
	case []float64:
		for _, f := range v {
			if math.IsNaN(f) {
				return true
			}
		}
	case []float32:
		for _, f := range v {
			if math.IsNaN(float64(f)) {
				return true
			}
		}
	case Matrix:
		return IsNan(v.DataP)
	case [4]Matrix:
		for n := 0; n < 4; n++ {
			if IsNan(v[n].DataP) {
				return true
			}
		}
	case [3]Matrix:
		for n := 0; n < 3; n++ {
			if IsNan(v[n].DataP) {
				return true
			}
		}
	case [2]Matrix:
		for n := 0; n < 2; n++ {
			if IsNan(v[n].DataP) {
				return true
			}
		}
	case Vector:
		return IsNan(v.DataP)
	}
	return false
}

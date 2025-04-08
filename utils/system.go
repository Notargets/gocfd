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

func IsNan(A any) {
	p := func(t bool) {
		if t {
			panic("NAN found")
		}
	}
	switch v := A.(type) {
	case float64:
		p(math.IsNaN(float64(v)))
	case float32:
		p(math.IsNaN(float64(v)))
	case []float64:
		for _, f := range v {
			p(math.IsNaN(f))
		}
	case []float32:
		for _, f := range v {
			p(math.IsNaN(float64(f)))
		}
	case Matrix:
		IsNan(v.DataP)
	case [4]Matrix:
		for n := 0; n < 4; n++ {
			IsNan(v[n].DataP)
		}
	case [3]Matrix:
		for n := 0; n < 3; n++ {
			IsNan(v[n].DataP)
		}
	case [2]Matrix:
		for n := 0; n < 2; n++ {
			IsNan(v[n].DataP)
		}
	case Vector:
		IsNan(v.DataP)
	}
}

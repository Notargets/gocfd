package utils

import (
	"fmt"
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

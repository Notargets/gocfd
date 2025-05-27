//go:build cgo
// +build cgo

package utils

/*
// Linux AMD64 - Full optimization with AVX support
#cgo linux,amd64 CFLAGS: -march=native -mavx -mavx2
#cgo linux,amd64 LDFLAGS: -lopenblas -llapacke -lgfortran -lm -lpthread

// Linux ARM64 - Native optimization without x86-specific flags
#cgo linux,arm64 CFLAGS: -march=native
#cgo linux,arm64 LDFLAGS: -lopenblas -llapacke -lgfortran -lm -lpthread

// Linux ARM (32-bit) - Conservative optimization
#cgo linux,arm CFLAGS: -march=native
#cgo linux,arm LDFLAGS: -lopenblas -llapacke -lgfortran -lm -lpthread

// macOS AMD64 - Full optimization with AVX support
#cgo darwin,amd64 CFLAGS: -march=native -mavx -mavx2
#cgo darwin,amd64 LDFLAGS: -lopenblas -llapacke -lgfortran -lm -lpthread

// macOS ARM64 (M1/M2/M3) - Native optimization for Apple Silicon
#cgo darwin,arm64 CFLAGS: -march=native
#cgo darwin,arm64 LDFLAGS: -lopenblas -llapacke -lgfortran -lm -lpthread

// Windows AMD64 - Conservative flags, may need different library names
#cgo windows,amd64 CFLAGS: -march=native
#cgo windows,amd64 LDFLAGS: -lm

// Windows ARM64 - Basic support
#cgo windows,arm64 CFLAGS: -march=native
#cgo windows,arm64 LDFLAGS: -lm

#include <cblas.h>
#include <lapacke.h>
*/
import "C"

import (
	"fmt"

	"gonum.org/v1/gonum/blas/blas64"
	netblas "gonum.org/v1/netlib/blas/netlib"
)

func init() {
	blas64.Use(netblas.Implementation{})
	fmt.Println("Using netlib to accelerate BLAS")
}

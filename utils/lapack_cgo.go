//go:build cgo
// +build cgo

package utils

/*
#cgo CFLAGS: -march=native -mavx -mavx2
#cgo LDFLAGS: -lopenblas -llapacke -lgfortran -lm -lpthread
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

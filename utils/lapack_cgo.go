//go:build cgo
// +build cgo

package utils

/*
#cgo CFLAGS: -O3 -march=native -mavx -mavx2 -mfma
#cgo LDFLAGS: -lopenblas -llapacke -lgfortran -lm -lpthread
#include <cblas.h>
#include <lapacke.h>
*/

import "C"

import (
	"gonum.org/v1/gonum/blas/blas64"
	"gonum.org/v1/gonum/lapack/lapack64"
	netblas "gonum.org/v1/netlib/blas/netlib"
	netlapack "gonum.org/v1/netlib/lapack/netlib"
)

func init() {
	blas64.Use(netblas.Implementation{})
	lapack64.Use(netlapack.Implementation{})
}

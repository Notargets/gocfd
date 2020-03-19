package utils

import (
	"gonum.org/v1/gonum/blas/blas64"
	"gonum.org/v1/gonum/mat"
)

// Row is a user-defined Row vector.
type Row []float64

// Dims, At and T minimally satisfy the mat.Matrix interface.
func (v Row) Dims() (r, c int)    { return 1, len(v) }
func (v Row) At(_, j int) float64 { return v[j] }
func (v Row) T() mat.Matrix       { return Column(v) }

// RawVector allows fast path computation with the vector.
func (v Row) RawVector() blas64.Vector {
	return blas64.Vector{N: len(v), Data: v, Inc: 1}
}

// Column is a user-defined Column vector.
type Column []float64

// Dims, At and T minimally satisfy the mat.Matrix interface.
func (v Column) Dims() (r, c int)    { return len(v), 1 }
func (v Column) At(i, _ int) float64 { return v[i] }
func (v Column) T() mat.Matrix       { return Row(v) }

// RawVector allows fast path computation with the vector.
func (v Column) RawVector() blas64.Vector {
	return blas64.Vector{N: len(v), Data: v, Inc: 1}
}

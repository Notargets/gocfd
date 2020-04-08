package utils

import (
	"gonum.org/v1/gonum/blas/blas64"
	"gonum.org/v1/gonum/mat"
)

type Vector struct {
	V *mat.VecDense
}

func NewVector(n int, dataO ...[]float64) Vector {
	var (
		data []float64
	)
	if len(dataO) != 0 {
		data = dataO[0]
	}
	return Vector{
		mat.NewVecDense(n, data),
	}
}

// Dims, At and T minimally satisfy the mat.Matrix interface.
func (v Vector) Dims() (r, c int)         { return v.V.Dims() }
func (v Vector) At(i, j int) float64      { return v.V.At(i, j) }
func (v Vector) T() mat.Matrix            { return v.V.T() }
func (v Vector) AtVec(i int) float64      { return v.V.AtVec(i) }
func (v Vector) RawVector() blas64.Vector { return v.V.RawVector() }
func (v Vector) SubVec(a, b Vector)       { v.V.SubVec(a.V, b.V) }
func (v Vector) Len() int                 { return v.V.Len() }

// Chainable (extended) methods
func (v Vector) Transpose() Matrix {
	var (
		nr, nc = v.V.Dims()
		m      *mat.Dense
	)
	m = mat.NewDense(nc, nr, v.V.RawVector().Data)
	return Matrix{m}
}

func (v Vector) Subtract(a Vector) Vector {
	var (
		data  = v.V.RawVector().Data
		dataA = a.V.RawVector().Data
	)
	for i := range data {
		data[i] -= dataA[i]
	}
	return v
}

func (v Vector) Subset(I Index) Vector {
	var (
		data  = v.V.RawVector().Data
		dataR = make([]float64, len(I))
		r     *mat.VecDense
	)
	for i, ind := range I {
		dataR[i] = data[ind]
	}
	r = mat.NewVecDense(len(dataR), dataR)
	return Vector{r}
}

func (v Vector) Scale(a float64) Vector {
	var (
		data = v.V.RawVector().Data
	)
	for i := range data {
		data[i] *= a
	}
	return v
}

func (v Vector) AddScalar(a float64) Vector {
	var (
		data = v.V.RawVector().Data
	)
	for i := range data {
		data[i] += a
	}
	return v
}

func (v Vector) Set(a float64) Vector {
	var (
		data = v.V.RawVector().Data
	)
	for i := range data {
		data[i] = a
	}
	return v
}

func (v Vector) Apply(f func(float64) float64) Vector {
	var (
		data = v.V.RawVector().Data
	)
	for i, val := range data {
		data[i] = f(val)
	}
	return v
}

func (v Vector) POW(p int) Vector {
	var (
		data = v.V.RawVector().Data
	)
	for i, val := range data {
		data[i] = POW(val, p)
	}
	return v
}

func (v Vector) Copy() Vector {
	return Vector{
		mat.VecDenseCopyOf(v.V),
	}
}

func (v Vector) Find(op EvalOp, target float64, abs bool) (r Vector) {
	r.V = VecFind(v.V, op, target, abs)
	return
}

// Non Chainable methods
func (v Vector) ToIndex() (I Index) {
	var (
		data = v.V.RawVector().Data
	)
	I = make(Index, v.Len())
	for i, val := range data {
		I[i] = int(val)
	}
	return
}

func (v Vector) Min() (min float64) {
	var (
		data = v.V.RawVector().Data
	)
	min = data[0]
	for _, val := range data {
		if val < min {
			min = val
		}
	}
	return
}

func (v Vector) Max() (max float64) {
	var (
		data = v.V.RawVector().Data
	)
	max = data[0]
	for _, val := range data {
		if val > max {
			max = val
		}
	}
	return
}

func (v Vector) ToMatrix() Matrix {
	var (
		m = mat.NewDense(v.V.Len(), 1, v.V.RawVector().Data)
	)
	return Matrix{m}
}

func (v Vector) Mul(w Vector) (A Matrix) {
	var (
		dataV = v.V.RawVector().Data
		dataW = w.V.RawVector().Data
		nr    = v.Len()
		nc    = w.Len()
		dataA []float64
	)
	A = NewMatrix(nr, nc)
	dataA = A.M.RawMatrix().Data
	for j := 0; j < nc; j++ {
		for i := 0; i < nr; i++ {
			ind := nc*i + j
			dataA[ind] = dataV[i] * dataW[j]
		}
	}
	return
}

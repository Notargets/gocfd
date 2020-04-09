package utils

import (
	"fmt"

	"gonum.org/v1/gonum/blas/blas64"

	"gonum.org/v1/gonum/mat"
)

type Matrix struct {
	M *mat.Dense
}

func NewMatrix(nr, nc int, dataO ...[]float64) Matrix {
	var (
		data []float64
	)
	if len(dataO) != 0 {
		data = dataO[0]
	}
	return Matrix{
		mat.NewDense(nr, nc, data),
	}
}

// Dims, At and T minimally satisfy the mat.Matrix interface.
func (m Matrix) Dims() (r, c int)          { return m.M.Dims() }
func (m Matrix) At(i, j int) float64       { return m.M.At(i, j) }
func (m Matrix) T() mat.Matrix             { return m.T() }
func (m Matrix) RawMatrix() blas64.General { return m.M.RawMatrix() }

// Chainable methods (extended)
func (m Matrix) Slice(I, K, J, L int) (R Matrix) {
	var (
		nrR   = K - I
		ncR   = L - J
		dataR = make([]float64, nrR*ncR)
		nr, _ = m.Dims()
		data  = m.M.RawMatrix().Data
	)
	for j := J; j < L; j++ {
		for i := I; i < K; i++ {
			ind := i + nr*j
			iR := i - I
			jR := j - J
			indR := iR + nrR*jR
			dataR[indR] = data[ind]
		}
	}
	R = NewMatrix(nrR, ncR, dataR)
	return
}

func (m Matrix) Copy() (R Matrix) {
	var (
		data   = m.M.RawMatrix().Data
		nr, nc = m.Dims()
		dataR  = make([]float64, nr*nc)
	)
	for i, val := range data {
		dataR[i] = val
	}
	R = NewMatrix(nr, nc, dataR)
	return
}

func (m Matrix) Set(i, j int, val float64) Matrix {
	m.M.Set(i, j, val)
	return m
}

func (m Matrix) Transpose() (R Matrix) {
	var (
		nr, nc = m.Dims()
		data   = m.M.RawMatrix().Data
	)
	R = NewMatrix(nc, nr)
	dataR := R.M.RawMatrix().Data
	for j := 0; j < nc; j++ {
		for i := 0; i < nr; i++ {
			ind := i + nr*j
			indR := i*nc + j
			dataR[indR] = data[ind]
		}
	}
	return
}

func (m Matrix) SetCol(j int, data []float64) Matrix {
	m.M.SetCol(j, data)
	return m
}

func (m Matrix) Mul(A Matrix) Matrix {
	var (
		nrM, ncM   = m.M.Dims()
		nrA, ncA   = A.M.Dims()
		r          = NewMatrix(nrM, ncA, nil)
		_, _, _, _ = nrM, ncM, nrA, ncA
	)
	r.M.Mul(m.M, A.M)
	return r
}

func (m Matrix) Add(A Matrix) Matrix {
	m.M.Add(m.M, A)
	return m
}

func (m Matrix) Subtract(a Matrix) Matrix {
	var (
		data  = m.M.RawMatrix().Data
		dataA = a.M.RawMatrix().Data
	)
	for i := range data {
		data[i] -= dataA[i]
	}
	return m
}

func (m Matrix) Subset(I Index) Matrix {
	/*
		Index should contain a list of indices into MI
		Note: native mat library matrix storage is in column traversal first (row-major) order
	*/
	var (
		Mr     = m.RawMatrix()
		nr, nc = m.Dims()
		data   = make([]float64, nr*nc)
		R      *mat.Dense
	)
	for _, ind := range I {
		data[ind] = Mr.Data[ind]
	}
	R = mat.NewDense(nr, nc, data)
	return Matrix{R}
}

func (m Matrix) Scale(a float64) Matrix {
	var (
		data = m.M.RawMatrix().Data
	)
	for i := range data {
		data[i] *= a
	}
	return m
}

func (m Matrix) AddScalar(a float64) Matrix {
	var (
		data = m.M.RawMatrix().Data
	)
	for i := range data {
		data[i] += a
	}
	return m
}

func (m Matrix) Apply(f func(float64) float64) Matrix {
	var (
		data = m.M.RawMatrix().Data
	)
	for i, val := range data {
		data[i] = f(val)
	}
	return m
}

func (m Matrix) POW(p int) Matrix {
	var (
		data = m.M.RawMatrix().Data
	)
	for i, val := range data {
		data[i] = POW(val, p)
	}
	return m
}

func (m Matrix) SliceRows(I Index) (R Matrix) {
	// RowIndices should contain a list of row indices into M
	var (
		nr, nc = m.Dims()
		nI     = len(I)
	)
	R = NewMatrix(nI, nc)
	for i, val := range I {
		if val > nr-1 || val < 0 {
			fmt.Printf("index out of bounds: index = %d, max_bounds = %d\n", val, nr-1)
			panic("unable to subset row from matrix")
		}
		R.M.SetRow(i, m.M.RawRowView(val))
	}
	return
}

func (m Matrix) ElementMultiply(A Matrix) Matrix {
	var (
		dataM = m.RawMatrix().Data
		dataA = A.RawMatrix().Data
	)
	for i, val := range dataA {
		dataM[i] *= val
	}
	return m
}

func (m Matrix) SubAssign(I Index, A Matrix) Matrix {
	var (
		dataM = m.RawMatrix().Data
		dataA = A.RawMatrix().Data
	)
	for _, ind := range I {
		dataM[ind] = dataA[ind]
	}
	return m
}

func (m Matrix) SubAssignScalar(I Index, val float64) Matrix {
	var (
		dataM = m.RawMatrix().Data
	)
	for _, ind := range I {
		dataM[ind] = val
	}
	return m
}

// Non chainable methods
func (m Matrix) IndexedAssign(I2 Index2D, Val Index) (err error) {
	var (
		data = m.RawMatrix().Data
	)
	if I2.Len != len(Val) {
		err = fmt.Errorf("length of index and values are not equal: len(I2) = %v, len(Val) = %v\n", I2.Len, len(Val))
		return
	}
	for i, val := range Val {
		data[I2.Ind[i]] = float64(val)
	}
	return
}

func (m Matrix) Inverse() (R Matrix, err error) {
	var (
		nr, nc = m.Dims()
	)
	R = NewMatrix(nr, nc)
	err = R.M.Inverse(m.M)
	return
}

func (m Matrix) Col(j int) Vector {
	var (
		data   = m.M.RawMatrix().Data
		nr, nc = m.M.Dims()
		vData  = make([]float64, nr)
	)
	for i := range vData {
		vData[i] = data[i*nc+j]
	}
	return Vector{
		mat.NewVecDense(nr, vData),
	}
}

func (m Matrix) Row(i int) Vector {
	var (
		data  = m.M.RawMatrix().Data
		_, nc = m.M.Dims()
		vData = make([]float64, nc)
	)
	for j := range vData {
		vData[j] = data[j+i*nc]
	}
	return Vector{
		mat.NewVecDense(nc, vData),
	}
}

func (m Matrix) Min() (min float64) {
	var (
		data = m.M.RawMatrix().Data
	)
	min = data[0]
	for _, val := range data {
		if val < min {
			min = val
		}
	}
	return
}

func (m Matrix) Max() (max float64) {
	var (
		data = m.M.RawMatrix().Data
	)
	max = data[0]
	for _, val := range data {
		if val > max {
			max = val
		}
	}
	return
}

func (m Matrix) Find(op EvalOp, val float64) (I Index2D) {
	var (
		nr, nc         = m.Dims()
		data           = m.RawMatrix().Data
		rowInd, colInd Index
	)
	switch op {
	case Equal:
		for j := 0; j < nc; j++ {
			for i := 0; i < nr; i++ {
				ind := i + nr*j
				if data[ind] == val {
					rowInd = append(rowInd, i)
					colInd = append(colInd, j)
				}
			}
		}
	case Less:
		for j := 0; j < nc; j++ {
			for i := 0; i < nr; i++ {
				ind := i + nr*j
				if data[ind] < val {
					rowInd = append(rowInd, i)
					colInd = append(colInd, j)
				}
			}
		}
	case LessOrEqual:
		for j := 0; j < nc; j++ {
			for i := 0; i < nr; i++ {
				ind := i + nr*j
				if data[ind] <= val {
					rowInd = append(rowInd, i)
					colInd = append(colInd, j)
				}
			}
		}
	case Greater:
		for j := 0; j < nc; j++ {
			for i := 0; i < nr; i++ {
				ind := i + nr*j
				if data[ind] > val {
					rowInd = append(rowInd, i)
					colInd = append(colInd, j)
				}
			}
		}
	case GreaterOrEqual:
		for j := 0; j < nc; j++ {
			for i := 0; i < nr; i++ {
				ind := i + nr*j
				if data[ind] >= val {
					rowInd = append(rowInd, i)
					colInd = append(colInd, j)
				}
			}
		}
	}
	I, _ = NewIndex2D(nr, nc, rowInd, colInd)
	return
}

func (m Matrix) SubsetVector(I Index) (V Vector) {
	var (
		Mr     = m.RawMatrix()
		nr, nc = m.Dims()
		data   = make([]float64, len(I))
	)
	for i, ind := range I {
		data[i] = Mr.Data[RowMajorToColMajor(nr, nc, ind)]
	}
	V = NewVector(len(I), data)
	return
}

func RowMajorToColMajor(nr, nc, ind int) (cind int) {
	// ind = i + nr * j
	// ind / nr = 0 + j
	j := ind / nr
	i := ind - nr*j
	cind = j + nc*i
	return
}

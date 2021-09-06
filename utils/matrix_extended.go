package utils

import (
	"fmt"
	"math"
	"sync"

	"gonum.org/v1/gonum/lapack/lapack64"

	"gonum.org/v1/gonum/blas/blas64"

	"gonum.org/v1/gonum/mat"
)

type Matrix struct {
	M        *mat.Dense
	readOnly bool
	name     string
	DataP    []float64
}

var (
	I, Zero Matrix
)

func init() {
	I = NewMatrix(1, 1, []float64{1.})
	Zero = NewMatrix(1, 1, []float64{0.})
}

func NewMatrix(nr, nc int, dataO ...[]float64) (R Matrix) {
	var (
		m        *mat.Dense
		dataArea []float64
	)
	if len(dataO) != 0 {
		if len(dataO[0]) < nr*nc {
			err := fmt.Errorf("mismatch in allocation: NewMatrix nr,nc = %v,%v, len(data[0]) = %v\n", nr, nc, len(dataO[0]))
			panic(err)
		}
		dataArea = dataO[0][0 : nr*nc]
	} else {
		dataArea = make([]float64, nr*nc)
	}
	m = mat.NewDense(nr, nc, dataArea)
	R = Matrix{
		M:        m,
		readOnly: false,
		name:     "unnamed - hint: pass a variable name to SetReadOnly()",
		DataP:    m.RawMatrix().Data,
	}
	return
}

// Dims, At and T minimally satisfy the mat.Matrix interface.
func (m Matrix) Dims() (r, c int)          { return m.M.Dims() }
func (m Matrix) At(i, j int) float64       { return m.M.At(i, j) }
func (m Matrix) T() mat.Matrix             { return m.T() }
func (m Matrix) RawMatrix() blas64.General { return m.M.RawMatrix() }
func (m Matrix) Data() []float64 {
	return m.RawMatrix().Data
}

// Utility
func (m *Matrix) IsEmpty() bool {
	if m.M == nil {
		return true
	}
	return false
}

// Chainable methods (extended)
func (m *Matrix) SetReadOnly(name ...string) Matrix {
	if len(name) != 0 {
		m.name = name[0]
	}
	m.readOnly = true
	return *m
}

func (m *Matrix) SetWritable() Matrix {
	m.readOnly = false
	return *m
}

func (m Matrix) Slice(I, K, J, L int) (R Matrix) { // Does not change receiver
	var (
		nrR   = K - I
		ncR   = L - J
		dataR = make([]float64, nrR*ncR)
		_, nc = m.Dims()
		data  = m.RawMatrix().Data
	)
	for j := J; j < L; j++ {
		for i := I; i < K; i++ {
			ind := i*nc + j
			iR := i - I
			jR := j - J
			indR := iR*ncR + jR
			dataR[indR] = data[ind]
		}
	}
	R = NewMatrix(nrR, ncR, dataR)
	return
}

func (m Matrix) Copy(RO ...Matrix) (R Matrix) { // Does not change receiver
	var (
		nr, nc = m.Dims()
	)
	R = getResultMatrix(nr, nc, RO)
	copy(R.DataP, m.DataP)
	return
}

func (m Matrix) Print(msgI ...string) (o string) {
	var (
		name = ""
	)
	if len(msgI) != 0 {
		name = msgI[0]
	}
	formatString := "%s = \n%8.5f\n"
	o = fmt.Sprintf(formatString, name, mat.Formatted(m.M, mat.Squeeze()))
	return
}

func (m Matrix) Transpose() (R Matrix) { // Does not change receiver
	var (
		nr, nc = m.Dims()
		data   = m.RawMatrix().Data
	)
	R = NewMatrix(nc, nr)
	dataR := R.M.RawMatrix().Data
	for j := 0; j < nc; j++ {
		for i := 0; i < nr; i++ {
			ind := i*nc + j
			indR := i + nr*j
			dataR[indR] = data[ind]
		}
	}
	return
}

func getResultMatrix(nr, nc int, RO []Matrix) (R Matrix) { // helper function - parses incoming optional destination argument
	if len(RO) != 0 {
		// Optional matrix for result is present
		R = RO[0]
		nrR, ncR := R.Dims()
		if nrR != nr || ncR != nc {
			panic("incorrect dimensions for provided result matrix")
		}
	} else {
		R = NewMatrix(nr, nc)
	}
	return
}

func (m Matrix) Mul(A Matrix, RO ...Matrix) (R Matrix) { // Does not change receiver
	var (
		nrM, ncM = m.Dims()
		nrA, ncA = A.M.Dims()
		nrR, ncR int
	)
	switch {
	case m.IsScalar():
		nrR, ncR = nrA, ncA
	case A.IsScalar():
		nrR, ncR = nrM, ncM
	default:
		nrR, ncR = nrM, ncA
	}
	R = getResultMatrix(nrR, ncR, RO)
	switch {
	case m.IsScalar():
		A.Copy(R)
		R.Scale(m.DataP[0])
	case A.IsScalar():
		m.Copy(R)
		R.Scale(A.DataP[0])
	default:
		R.M.Mul(m.M, A.M)
	}
	return R
}
func (m Matrix) MulParallel(A Matrix, nP int) (R Matrix) { // Does not change receiver
	var (
		nrM, _   = m.Dims()
		nrA, ncA = A.M.Dims()
		wg       = sync.WaitGroup{}
		aD       = A.DataP
	)
	MulExt := func(m, A Matrix, dataO ...[]float64) (R Matrix) {
		var (
			nrM, _ = m.Dims()
			_, ncA = A.M.Dims()
		)
		if len(dataO) != 0 {
			R = NewMatrix(nrM, ncA, dataO[0])
		} else {
			R = NewMatrix(nrM, ncA)
		}
		R.M.Mul(m.M, A.M)
		return R
	}
	if nP > ncA {
		nP = ncA
	}
	R = NewMatrix(nrM, ncA)
	rD := R.DataP
	ncAChunk := Split1DMaxChunk(ncA, nP)
	subAChunkSize := nrA * ncAChunk
	subRChunkSize := nrM * ncAChunk
	subAstorage := make([]float64, nP*subAChunkSize)
	subRstorage := make([]float64, nP*subRChunkSize)
	for n := 0; n < nP; n++ {
		ind, end := Split1D(ncA, nP, n)
		ncSubA := end - ind
		wg.Add(1)
		go func(ind, end, ncSubA, n int) {
			subA := NewMatrix(nrA, ncSubA, subAstorage[n*subAChunkSize:])
			sAD := subA.DataP
			for j := 0; j < nrA; j++ {
				var ii int
				for i := ind; i < end; i++ {
					sAD[ii+ncSubA*j] = aD[i+ncA*j]
					ii++
				}
			}
			subR := MulExt(m, subA, subRstorage[n*subRChunkSize:])
			sRD := subR.DataP
			for j := 0; j < nrM; j++ {
				var ii int
				for i := ind; i < end; i++ {
					rD[i+ncA*j] = sRD[ii+ncSubA*j]
					ii++
				}
			}
			wg.Done()
		}(ind, end, ncSubA, n)
	}
	wg.Wait()
	return R
}

func getDimensions(m, A Matrix) (nr, nc int) {
	var (
		nrA, ncA = A.Dims()
	)
	nr, nc = m.Dims()
	switch {
	case m.IsScalar():
		nr, nc = A.Dims() // upscale m (scalar) to A
	case A.IsScalar():
		nrA, ncA = m.Dims() // upscale A (scalar) to m
	}
	if nrA != nr || ncA != nc {
		panic("dimensions of matrices do not match")
	}
	return
}

func (m *Matrix) upscale(nr, nc int) {
	if m.IsScalar() {
		val := m.DataP[0]
		RR := NewMatrix(nr, nc)
		for i := 0; i < nr; i++ {
			RR.Set(i, i, val)
		}
		m.M, m.DataP = RR.M, RR.DataP
	}
}

func add(m, A Matrix, subtract bool, RO []Matrix) (R Matrix) { // Changes receiver, optionally does not
	var (
		nr, nc = getDimensions(m, A)
		mult   = 1.
	)
	if subtract {
		mult = -1.
	}
	if len(RO) != 0 {
		R = getResultMatrix(nr, nc, RO)
		if m.IsScalar() {
			val := m.DataP[0]
			for i := 0; i < nr; i++ {
				R.Set(i, i, val)
			}
		} else {
			for i, val := range m.DataP {
				R.DataP[i] = val
			}
		}
	} else {
		m.checkWritable()
		m.upscale(nr, nc) // Must upscale m to match A
		R = m
	}
	switch {
	case A.IsScalar():
		AVal := mult * A.DataP[0]
		for i := 0; i < nr; i++ {
			R.Set(i, i, R.At(i, i)+AVal)
		}
	default:
		for i := range A.DataP {
			R.DataP[i] += mult * A.DataP[i]
		}
	}
	return
}

func (m Matrix) Add(A Matrix, RO ...Matrix) Matrix { // Changes receiver optionally
	return add(m, A, false, RO)
}

func (m Matrix) Subtract(A Matrix, RO ...Matrix) Matrix { // Changes receiver optionally
	return add(m, A, true, RO)
}

func (m Matrix) AddScalar(a float64) Matrix { // Changes receiver
	var (
		data = m.RawMatrix().Data
	)
	m.checkWritable()
	for i := range data {
		data[i] += a
	}
	return m
}

func (m Matrix) Subset(I Index, newDims ...int) (R Matrix) { // Does not change receiver
	/*
		Index should contain a list of indices into MI
		Note: native mat library matrix storage is in column traversal first (row-major) order
	*/
	var (
		Mr           = m.RawMatrix()
		nr, nc       = m.Dims()
		data         []float64
		nrNew, ncNew = nr, nc
	)
	if len(newDims) != 0 {
		nrNew, ncNew = newDims[0], newDims[1]
	}
	data = make([]float64, nrNew*ncNew)
	for i, ind := range I {
		data[i] = Mr.Data[ind]
	}
	return NewMatrix(nrNew, ncNew, data)
}

func (m Matrix) SliceRows(I Index) (R Matrix) { // Does not change receiver
	// I should contain a list of row indices into M
	var (
		nr, nc   = m.Dims()
		nI       = len(I)
		maxIndex = nr - 1
	)
	R = NewMatrix(nI, nc)
	for iNewRow, i := range I {
		if i > maxIndex || i < 0 {
			fmt.Printf("index out of bounds: index = %d, max_bounds = %d\n", i, maxIndex)
			panic("unable to subset rows from matrix")
		}
		R.M.SetRow(iNewRow, m.M.RawRowView(i))
	}
	return
}

func (m Matrix) SliceCols(I Index) (R Matrix) { // Does not change receiver
	// I should contain a list of column indices into M
	var (
		nr, nc   = m.Dims()
		maxIndex = nc - 1
		nI       = len(I)
		dataM    = m.RawMatrix().Data
		colData  = make([]float64, nr)
	)
	R = NewMatrix(nr, nI)
	for jNewCol, j := range I {
		if j > maxIndex || j < 0 {
			fmt.Printf("index out of bounds: index = %d, max_bounds = %d\n", j, maxIndex)
			panic("unable to subset columns from matrix")
		}
		var ind int
		for i := 0; i < nr; i++ {
			ind = i*nc + j
			colData[i] = dataM[ind]
		}
		R.M.SetCol(jNewCol, colData)
	}
	return
}

func (m Matrix) Set(i, j int, val float64) Matrix { // Changes receiver
	var (
		nr, nc = m.Dims()
	)
	i, j = lim(i, nr), lim(j, nc)
	m.checkWritable()
	m.M.Set(i, j, val)
	return m
}

func (m Matrix) SetRange(i1, i2, j1, j2 int, val float64) Matrix { // Changes receiver
	var (
		nr, nc = m.Dims()
		data   = m.RawMatrix().Data
	)
	m.checkWritable()
	i1, i2, j1, j2 = limRange(i1, i2, j1, j2, nr, nc)
	for i := i1; i < i2; i++ {
		for j := j1; j < j2; j++ {
			ind := i*nc + j
			data[ind] = val
		}
	}
	return m
}

func (m Matrix) SetCol(j int, data []float64) Matrix { // Changes receiver
	var (
		_, nc = m.Dims()
	)
	j = lim(j, nc)
	m.checkWritable()
	m.M.SetCol(j, data)
	return m
}

func (m Matrix) SetColFrom(col, rowFrom int, data []float64) Matrix { // Changes receiver
	var (
		nr, nc = m.Dims()
	)
	col = lim(col, nc)
	m.checkWritable()
	if len(data)+rowFrom > nr {
		panic(
			fmt.Errorf("row length exceeded, max is %d, have %d", nr, rowFrom+len(data)))
	}
	for i, val := range data {
		m.M.Set(i+rowFrom, col, val)
	}
	return m
}

func (m Matrix) AssignColumns(I Index, A Matrix) Matrix { // Changes receiver
	var (
		_, nc = m.Dims()
	)
	// Assigns columns in M to columns indexed sequentially from Ainv
	m.checkWritable()
	for i, j := range I {
		if j > nc-1 {
			err := fmt.Errorf("bad index value: %v exceeds bounds of %v", j, nc)
			panic(err)
		}
		m.SetCol(j, A.Col(i).RawVector().Data)
	}
	return m
}

func (m Matrix) Assign(I Index, AI interface{}) Matrix {
	if err := m.IndexedAssign(I, AI); err != nil {
		panic(err)
	}
	return m
}

func (m Matrix) Range(RangeO ...interface{}) (r []float64) {
	var (
		nr, nc = m.Dims()
		I      Index
		data   = m.Data()
	)
	I = expandRangeO(nr, nc, RangeO)
	r = make([]float64, len(I))
	for i, ind := range I {
		r[i] = data[ind]
	}
	return
}

func (m Matrix) Equate(ValuesI interface{}, RangeO ...interface{}) {
	var (
		nr, nc = m.Dims()
		I      Index
		Values []float64
		nVal   int
	)
	I = expandRangeO(nr, nc, RangeO)
	nVal = len(I)
	Values = expandValues(nVal, ValuesI)
	m.Assign(I, Values)
}

func expandValues(nVal int, ValuesI interface{}) (vals []float64) {
	switch values := ValuesI.(type) {
	case []float64:
		vals = values
		if len(vals) != nVal {
			goto FAIL
		}
	case []int:
		vals = make([]float64, nVal)
		for i := range vals {
			vals[i] = float64(values[i])
		}
	case Index:
		vals = make([]float64, nVal)
		for i := range vals {
			vals[i] = float64(values[i])
		}
	case int:
		vals = make([]float64, nVal)
		for i := range vals {
			vals[i] = float64(values)
		}
	case float64:
		vals = make([]float64, nVal)
		for i := range vals {
			vals[i] = values
		}
	case Matrix:
		vals = values.Data()
		if len(vals) != nVal {
			goto FAIL
		}
	case Vector:
		vals = values.Data()
		if len(vals) != nVal {
			goto FAIL
		}
	}
	return
FAIL:
	panic(fmt.Errorf("number of values not equal to index"))
}

func expandRangeO(nr, nc int, RangeO []interface{}) (I Index) {
	var (
		err error
		I2D Index2D
	)
	switch len(RangeO) {
	case 1:
		I = expandRangeI(nr, RangeO[0])
	case 2:
		I1 := expandRangeI(nr, RangeO[0])
		I2 := expandRangeI(nc, RangeO[1])
		if I2D, err = NewIndex2D(nr, nc, I1, I2, true); err != nil {
			panic(err)
		}
		I = I2D.ToIndex()
	default:
		panic(fmt.Errorf("only 1D and 2D ranges supported"))
	}
	return
}

func expandRangeI(max int, RangeI interface{}) (I Index) {
	switch val := RangeI.(type) {
	case []int:
		I = val
	case Index:
		I = val
	case []float64:
		I = make(Index, len(val))
		for i, val := range val {
			I[i] = int(val)
		}
	case int:
		I = make(Index, 1)
		I[0] = val
	case string:
		r1 := NewR1(max)
		I = r1.Range(val)
	case Vector:
		I = expandRangeI(max, val.Data())
	case Matrix:
		I = expandRangeI(max, val.Data())
	}
	for _, val := range I {
		if val > max {
			panic(fmt.Errorf("max value %d exceeded by index value %d", max, val))
		}
	}
	return
}

//func (m Matrix) AssignVector(I Index, Ainv Vector) Matrix { // Changes receiver
func (m Matrix) AssignVector(I Index, AI interface{}) Matrix { // Changes receiver
	// Assigns values indexed into M using values sequentially from Vector Ainv
	var (
		dataM = m.RawMatrix().Data
	)
	m.checkWritable()
	switch A := AI.(type) {
	case Vector:
		dataA := A.RawVector().Data
		for i, ind := range I {
			dataM[ind] = dataA[i]
		}
	case Matrix:
		dataA := A.RawMatrix().Data
		for i, ind := range I {
			dataM[ind] = dataA[i]
		}
	case Index:
		for i, ind := range I {
			dataM[ind] = float64(A[i])
		}
	}
	return m
}

func (m Matrix) Scale(a float64) Matrix { // Changes receiver
	var (
		data = m.RawMatrix().Data
	)
	m.checkWritable()
	for i := range data {
		data[i] *= a
	}
	return m
}

func (m Matrix) Apply(f func(float64) float64) Matrix { // Changes receiver
	var (
		data = m.RawMatrix().Data
	)
	m.checkWritable()
	for i, val := range data {
		data[i] = f(val)
	}
	return m
}

func (m Matrix) Apply2(A Matrix, f func(float64, float64) float64) Matrix { // Changes receiver
	var (
		dataM = m.RawMatrix().Data
		dataA = A.RawMatrix().Data
	)
	m.checkWritable()
	for i, val := range dataM {
		dataM[i] = f(val, dataA[i])
	}
	return m
}

func (m Matrix) Apply3(A, B Matrix, f func(float64, float64, float64) float64) Matrix { // Changes receiver
	var (
		dataM  = m.RawMatrix().Data
		dA, dB = A.RawMatrix().Data, B.RawMatrix().Data
	)
	m.checkWritable()
	for i, val := range dataM {
		dataM[i] = f(val, dA[i], dB[i])
	}
	return m
}

func (m Matrix) Apply3Parallel(A, B Matrix, f func(float64, float64, float64) float64, nP int) Matrix { // Changes receiver
	var (
		dataM  = m.RawMatrix().Data
		dA, dB = A.RawMatrix().Data, B.RawMatrix().Data
		wg     = sync.WaitGroup{}
		l      = len(dataM)
	)
	m.checkWritable()
	for n := 0; n < nP; n++ {
		ind, end := Split1D(l, nP, n)
		wg.Add(1)
		go func(ind, end int) {
			for i := ind; i < end; i++ {
				val := dataM[i]
				dataM[i] = f(val, dA[i], dB[i])
			}
			wg.Done()
		}(ind, end)
	}
	wg.Wait()
	return m
}

func (m Matrix) Apply4(A, B, C Matrix, f func(float64, float64, float64, float64) float64) Matrix { // Changes receiver
	var (
		dataM      = m.RawMatrix().Data
		dA, dB, dC = A.RawMatrix().Data, B.RawMatrix().Data, C.RawMatrix().Data
	)
	m.checkWritable()
	for i, val := range dataM {
		dataM[i] = f(val, dA[i], dB[i], dC[i])
	}
	return m
}

func (m Matrix) Apply4Parallel(A, B, C Matrix, f func(float64, float64, float64, float64) float64, nP int) Matrix { // Changes receiver
	var (
		dataM      = m.RawMatrix().Data
		dA, dB, dC = A.RawMatrix().Data, B.RawMatrix().Data, C.RawMatrix().Data
		wg         = sync.WaitGroup{}
		l          = len(dataM)
	)
	m.checkWritable()
	for n := 0; n < nP; n++ {
		ind, end := Split1D(l, nP, n)
		wg.Add(1)
		go func(ind, end int) {
			for i := ind; i < end; i++ {
				val := dataM[i]
				dataM[i] = f(val, dA[i], dB[i], dC[i])
			}
			wg.Done()
		}(ind, end)
	}
	wg.Wait()
	return m
}

func (m Matrix) Apply5Parallel(A, B, C, D Matrix, f func(float64, float64, float64, float64, float64) float64, nP int) Matrix { // Changes receiver
	var (
		dataM          = m.RawMatrix().Data
		dA, dB, dC, dD = A.RawMatrix().Data, B.RawMatrix().Data, C.RawMatrix().Data, D.RawMatrix().Data
		wg             = sync.WaitGroup{}
		l              = len(dataM)
	)
	m.checkWritable()
	for n := 0; n < nP; n++ {
		ind, end := Split1D(l, nP, n)
		wg.Add(1)
		go func(ind, end int) {
			for i := ind; i < end; i++ {
				val := dataM[i]
				dataM[i] = f(val, dA[i], dB[i], dC[i], dD[i])
			}
			wg.Done()
		}(ind, end)
	}
	wg.Wait()
	return m
}

func (m Matrix) Apply8(A, B, C, D, E, F, G Matrix, f func(float64, float64, float64, float64, float64, float64, float64, float64) float64) Matrix { // Changes receiver
	var (
		dataM                      = m.RawMatrix().Data
		dA, dB, dC, dD, dE, dF, dG = A.RawMatrix().Data, B.RawMatrix().Data, C.RawMatrix().Data, D.RawMatrix().Data, E.RawMatrix().Data, F.RawMatrix().Data, G.RawMatrix().Data
	)
	m.checkWritable()
	for i, val := range dataM {
		dataM[i] = f(val, dA[i], dB[i], dC[i], dD[i], dE[i], dF[i], dG[i])
	}
	return m
}

func (m Matrix) POW(p int) Matrix { // Changes receiver
	var (
		data = m.RawMatrix().Data
	)
	m.checkWritable()
	for i, val := range data {
		data[i] = POW(val, p)
	}
	return m
}

func (m Matrix) ElMul(A Matrix) Matrix { // Changes receiver
	var (
		dataM = m.RawMatrix().Data
		dataA = A.RawMatrix().Data
	)
	m.checkWritable()
	for i, val := range dataA {
		dataM[i] *= val
	}
	return m
}

func (m Matrix) ElDiv(A Matrix) Matrix { // Changes receiver
	var (
		dataM = m.RawMatrix().Data
		dataA = A.RawMatrix().Data
	)
	m.checkWritable()
	for i, val := range dataA {
		dataM[i] /= val
	}
	return m
}

func (m Matrix) AssignScalar(I Index, val float64) Matrix { // Changes receiver
	var (
		dataM = m.RawMatrix().Data
	)
	m.checkWritable()
	for _, ind := range I {
		dataM[ind] = val
	}
	return m
}

func (m Matrix) LUSolve(b Matrix) (X Matrix) {
	var (
		err error
	)
	lu := &mat.LU{}
	lu.Factorize(m)
	X = NewMatrix(b.Dims())
	if err = lu.SolveTo(X.M, false, b.M); err != nil {
		panic(err)
	}
	return
}

// Non chainable methods
func (m Matrix) IsScalar() bool {
	var (
		nr, nc = m.Dims()
	)
	return nr*nc == 1
}

func (m Matrix) IndexedAssign(I Index, ValI interface{}) (err error) { // Changes receiver
	return IndexedAssign(m, I, ValI)
}

func (m Matrix) Inverse() (R Matrix, err error) {
	var (
		nr, nc      = m.Dims()
		errSingular = fmt.Errorf("unable to invert, matrix is singular")
	)
	R = m.Copy()
	if m.IsScalar() {
		if m.DataP[0] == 0. {
			err = errSingular
		} else {
			R.DataP[0] = 1. / m.DataP[0]
		}
		return
	}
	iPiv := make([]int, nr)
	if ok := lapack64.Getrf(R.RawMatrix(), iPiv); !ok {
		err = errSingular
		return
	}
	work := make([]float64, nr*nc)
	if ok := lapack64.Getri(R.RawMatrix(), iPiv, work, nr*nc); !ok {
		err = errSingular
	}
	return
}

func (m Matrix) Col(j int) Vector {
	var (
		data   = m.RawMatrix().Data
		nr, nc = m.Dims()
		vData  = make([]float64, nr)
	)
	j = lim(j, nc)
	for i := range vData {
		vData[i] = data[i*nc+j]
	}
	return NewVector(nr, vData)
}

func (m Matrix) Row(i int) Vector {
	var (
		data   = m.RawMatrix().Data
		nr, nc = m.Dims()
		vData  = make([]float64, nc)
	)
	i = lim(i, nr)
	for j := range vData {
		vData[j] = data[j+i*nc]
	}
	return NewVector(nc, vData)
}

func (m Matrix) Min() (min float64) {
	var (
		data = m.RawMatrix().Data
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
		data = m.RawMatrix().Data
	)
	max = data[0]
	for _, val := range data {
		if val > max {
			max = val
		}
	}
	return
}

func (m Matrix) Find(op EvalOp, val float64, abs bool) (I Index) {
	var (
		nr, nc = m.Dims()
		data   = m.RawMatrix().Data
	)
	var target float64
	switch op {
	case Equal:
		for j := 0; j < nc; j++ {
			for i := 0; i < nr; i++ {
				ind := i*nc + j
				target = data[ind]
				if abs {
					target = math.Abs(target)
				}
				if target == val {
					I = append(I, ind)
				}
			}
		}
	case Less:
		for j := 0; j < nc; j++ {
			for i := 0; i < nr; i++ {
				ind := i*nc + j
				target = data[ind]
				if abs {
					target = math.Abs(target)
				}
				if target < val {
					I = append(I, ind)
				}
			}
		}
	case LessOrEqual:
		for j := 0; j < nc; j++ {
			for i := 0; i < nr; i++ {
				ind := i*nc + j
				target = data[ind]
				if abs {
					target = math.Abs(target)
				}
				if target <= val {
					I = append(I, ind)
				}
			}
		}
	case Greater:
		for j := 0; j < nc; j++ {
			for i := 0; i < nr; i++ {
				ind := i*nc + j
				target = data[ind]
				if abs {
					target = math.Abs(target)
				}
				if target > val {
					I = append(I, ind)
				}
			}
		}
	case GreaterOrEqual:
		for j := 0; j < nc; j++ {
			for i := 0; i < nr; i++ {
				ind := i*nc + j
				target = data[ind]
				if abs {
					target = math.Abs(target)
				}
				if target >= val {
					I = append(I, ind)
				}
			}
		}
	}
	return
}

func (m Matrix) SubsetVector(I Index) (V Vector) {
	var (
		Mr = m.RawMatrix()
		//nr, nc = m.Dims()
		data = make([]float64, len(I))
	)
	for i, ind := range I {
		data[i] = Mr.Data[ind]
	}
	V = NewVector(len(I), data)
	return
}

func (m Matrix) SumRows() (V Vector) {
	/*
		Calculates the sum of each row to form the output vector, one result per row
	*/
	var (
		nr, nc = m.Dims()
		dataM  = m.RawMatrix().Data
		dataV  = make([]float64, nr)
	)
	for i := 0; i < nr; i++ {
		for j := 0; j < nc; j++ {
			dataV[i] += dataM[i*nc+j]
		}
	}
	V = NewVector(nr, dataV)
	return
}

func (m Matrix) SumCols() (V Vector) {
	/*
		Calculates the sum of each column to form the output vector, one result per column
	*/
	var (
		nr, nc = m.Dims()
		dataM  = m.RawMatrix().Data
		dataV  = make([]float64, nc)
	)
	for i := 0; i < nr; i++ {
		for j := 0; j < nc; j++ {
			dataV[j] += dataM[i*nc+j]
		}
	}
	V = NewVector(nc, dataV)
	return
}

func (m Matrix) checkWritable() {
	if m.readOnly {
		err := fmt.Errorf("attempt to write to a read only matrix named: \"%v\"", m.name)
		panic(err)
	}
}

func RowMajorToColMajor(nr, nc, ind int) (cind int) {
	// ind = i + nr * j
	// ind / nr = 0 + j
	j := ind / nr
	i := ind - nr*j
	cind = j + nc*i
	return
}

func lim(i, imax int) int {
	if i < 0 {
		return imax + i // Support indexing from end, -1 is imax
	}
	return i
}

func limLoop(ib, ie, imax int) (ibeg, iend int) {
	if ib < 0 {
		ibeg = imax + ib
	} else {
		ibeg = ib
	}
	if ie < 0 {
		iend = imax + ie + 1 // Support indexing from end, -1 is imax
	} else {
		iend = ie + 1
	}
	return
}

func limRange(i1, i2, j1, j2, nr, nc int) (ii1, ii2, jj1, jj2 int) {
	ii1, ii2 = limLoop(i1, i2, nr)
	jj1, jj2 = limLoop(j1, j2, nc)
	return ii1, ii2, jj1, jj2
}

func IndexedAssign(mI interface{}, I Index, ValI interface{}) (err error) { // Changes receiver
	var (
		temp []float64
	)
	switch Val := ValI.(type) {
	case []float64:
		temp = Val
	case Matrix:
		temp = Val.DataP
	case Index:
		temp = make([]float64, len(I))
		for i, val := range Val {
			temp[i] = float64(val)
		}
	}
	if len(I) != len(temp) {
		err = fmt.Errorf("length of index and values are not equal: len(I) = %v, len(Val) = %v\n", len(I), len(temp))
		return
	}
	switch m := mI.(type) {
	case Matrix:
		var data = m.RawMatrix().Data
		for i, val := range temp {
			data[I[i]] = val
		}
	case DOK:
		//_, nc := m.Dims()
		nr, _ := m.Dims()
		for ii, val := range temp {
			// DOK is stored column major, while the composed Index for the range is row-major, so we convert it
			i, j := indexToIJColMajor(I[ii], nr)
			m.M.Set(i, j, val)
		}
	case CSR:
		nr, _ := m.Dims()
		for ii, val := range temp {
			// CSR is stored column major, while the composed Index for the range is row-major, so we convert it
			i, j := indexToIJColMajor(I[ii], nr)
			m.M.Set(i, j, val)
		}
	}
	return
}

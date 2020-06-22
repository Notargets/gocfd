package utils

import (
	"fmt"

	"github.com/james-bowman/sparse"
	"github.com/james-bowman/sparse/blas"
	"gonum.org/v1/gonum/mat"
)

type DOK struct {
	M        *sparse.DOK
	readOnly bool
	name     string
}

func NewDOK(nr, nc int) (R DOK) {
	R = DOK{
		sparse.NewDOK(nr, nc),
		false,
		"unnamed - hint: pass a variable name to SetReadOnly()",
	}
	return
}

// Dims, At and T minimally satisfy the mat.Matrix interface.
func (m DOK) Dims() (r, c int)              { return m.M.Dims() }
func (m DOK) At(i, j int) float64           { return m.M.At(i, j) }
func (m DOK) T() mat.Matrix                 { return m.T() }
func (m DOK) RawMatrix() *blas.SparseMatrix { return m.M.RawMatrix() }
func (m DOK) Data() []float64 {
	return m.RawMatrix().Data
}

func (m DOK) Assign(I Index, AI interface{}) DOK {
	if err := m.IndexedAssign(I, AI); err != nil {
		panic(err)
	}
	return m
}

func (m DOK) Equate(ValuesI interface{}, RangeO ...interface{}) {
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

func (m DOK) IndexedAssign(I Index, ValI interface{}) (err error) { // Changes receiver
	var (
		temp  []float64
		_, nc = m.Dims()
	)
	m.checkWritable()
	switch Val := ValI.(type) {
	case []float64:
		temp = Val
	case Matrix:
		temp = Val.Data()
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
	for ii, val := range temp {
		// DOK is stored column major, while the composed Index for the range is row-major, so we convert it
		i, j := indexToIJColMajor(I[ii], nc)
		m.M.Set(i, j, val)
	}
	return
}

func (m DOK) checkWritable() {
	if m.readOnly {
		err := fmt.Errorf("attempt to write to a read only matrix named: \"%v\"", m.name)
		panic(err)
	}
}

func (m DOK) ToCSR() CSR {
	return CSR{
		M:        m.M.ToCSR(),
		readOnly: m.readOnly,
		name:     m.name,
	}
}

type CSR struct {
	M        *sparse.CSR
	readOnly bool
	name     string
}

func NewCSR(nr, nc int) (R CSR) {
	R = CSR{
		// TODO: This creates an unusable CSR, as the newly created matrix has no storage
		sparse.NewCSR(nr, nc, nil, nil, nil),
		false,
		"unnamed - hint: pass a variable name to SetReadOnly()",
	}
	return
}

// Dims, At and T minimally satisfy the mat.Matrix interface.
func (m CSR) Dims() (r, c int)              { return m.M.Dims() }
func (m CSR) At(i, j int) float64           { return m.M.At(i, j) }
func (m CSR) T() mat.Matrix                 { return m.T() }
func (m CSR) RawMatrix() *blas.SparseMatrix { return m.M.RawMatrix() }
func (m CSR) Data() []float64 {
	return m.RawMatrix().Data
}

func (m CSR) Assign(I Index, AI interface{}) CSR {
	if err := m.IndexedAssign(I, AI); err != nil {
		panic(err)
	}
	return m
}

func (m CSR) Equate(ValuesI interface{}, RangeO ...interface{}) {
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

func (m CSR) IndexedAssign(I Index, ValI interface{}) (err error) { // Changes receiver
	var (
		temp  []float64
		_, nc = m.Dims()
	)
	m.checkWritable()
	switch Val := ValI.(type) {
	case []float64:
		temp = Val
	case Matrix:
		temp = Val.Data()
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
	for ii, val := range temp {
		// CSR is stored column major, while the composed Index for the range is row-major, so we convert it
		i, j := indexToIJColMajor(I[ii], nc)
		fmt.Println("i, j = ", i, j)
		m.M.Set(i, j, val)
	}
	return
}

func (m CSR) checkWritable() {
	if m.readOnly {
		err := fmt.Errorf("attempt to write to a read only matrix named: \"%v\"", m.name)
		panic(err)
	}
}

func indexToIJ(ind, nc int) (i, j int) {
	//ind = j + nc*(i)
	i = ind / nc
	j = ind - i*nc
	return
}

func indexToIJColMajor(ind, nr int) (i, j int) {
	//ind = i + nr*(j)
	j = ind / nr
	i = ind - j*nr
	return
}

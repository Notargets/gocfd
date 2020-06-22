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

func (m DOK) checkWritable() {
	if m.readOnly {
		err := fmt.Errorf("attempt to write to a read only matrix named: \"%v\"", m.name)
		panic(err)
	}
}

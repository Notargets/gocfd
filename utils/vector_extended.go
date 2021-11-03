package utils

import (
	"fmt"
	"math"

	"gonum.org/v1/gonum/blas/blas64"
	"gonum.org/v1/gonum/mat"
)

type Vector struct {
	V     *mat.VecDense
	DataP []float64
}

func NewVector(n int, dataO ...[]float64) Vector {
	var (
		data []float64
	)
	if len(dataO) != 0 {
		data = dataO[0]
	}
	v := mat.NewVecDense(n, data)
	return Vector{
		V:     v,
		DataP: v.RawVector().Data,
	}
}

func NewVectorConstant(n int, val float64) Vector {
	var (
		dataV = ConstArray(n, val)
	)
	return NewVector(n, dataV)
}

// Dims, At and T minimally satisfy the mat.Matrix interface.
func (v Vector) Dims() (r, c int)         { return v.V.Dims() }
func (v Vector) At(i, j int) float64      { return v.V.At(i, j) }
func (v Vector) T() mat.Matrix            { return v.V.T() }
func (v Vector) AtVec(i int) float64      { return v.V.AtVec(i) }
func (v Vector) RawVector() blas64.Vector { return v.V.RawVector() }
func (v Vector) SubVec(a, b Vector)       { v.V.SubVec(a.V, b.V) }
func (v Vector) Len() int {
	if v.V != nil {
		return v.V.Len()
	} else {
		return 0
	}
}
func (v Vector) Data() []float64 {
	return v.DataP
}

// Chainable (extended) methods
func (v Vector) Linspace(begin, end float64) Vector {
	var (
		data = v.Data()
		rge  = (end - begin) / float64(v.Len()-1)
	)
	for i := range data {
		data[i] = begin + float64(i)*rge
	}
	return v
}

func (v Vector) Subtract(a Vector) Vector {
	var (
		data  = v.RawVector().Data
		dataA = a.RawVector().Data
	)
	for i := range data {
		data[i] -= dataA[i]
	}
	return v
}

func (v Vector) ElMul(a Vector) Vector {
	var (
		data  = v.RawVector().Data
		dataA = a.RawVector().Data
	)
	for i := range data {
		data[i] *= dataA[i]
	}
	return v
}

func (v Vector) Zip(op EvalOp, abs bool, A Vector) (R Vector) {
	/*
		Choose each element from v and a based on:
		if v op a == true {choose v} else {choose a}
	*/
	var (
		N  = v.Len()
		vD = v.RawVector().Data
		aD = A.RawVector().Data
		rD = make([]float64, N)
	)
	if N != A.Len() {
		err := fmt.Errorf("dimension mismatch: Zip receiver: %v Ainv: %v\n", N, A.Len())
		panic(err)
	}
	var target float64
	switch op {
	case Equal:
		for i, val := range vD {
			target = aD[i]
			if abs {
				val = math.Abs(val)
				target = math.Abs(target)
			}
			if val == aD[i] {
				rD[i] = val
			} else {
				rD[i] = target
			}
		}
	case Less:
		for i, val := range vD {
			target = aD[i]
			if abs {
				val = math.Abs(val)
				target = math.Abs(target)
			}
			if val < target {
				rD[i] = val
			} else {
				rD[i] = target
			}
		}
	case Greater:
		for i, val := range vD {
			target = aD[i]
			if abs {
				val = math.Abs(val)
				target = math.Abs(target)
			}
			if val > target {
				rD[i] = val
			} else {
				rD[i] = target
			}
		}
	case LessOrEqual:
		for i, val := range vD {
			target = aD[i]
			if abs {
				val = math.Abs(val)
				target = math.Abs(target)
			}
			if val <= target {
				rD[i] = val
			} else {
				rD[i] = target
			}
		}
	case GreaterOrEqual:
		for i, val := range vD {
			target = aD[i]
			if abs {
				val = math.Abs(val)
				target = math.Abs(target)
			}
			if val >= target {
				rD[i] = val
			} else {
				rD[i] = target
			}
		}
	}
	R = NewVector(N, rD)
	return R
}

func (v Vector) ElDiv(a Vector) Vector {
	var (
		data  = v.RawVector().Data
		dataA = a.RawVector().Data
	)
	for i := range data {
		data[i] /= dataA[i]
	}
	return v
}

func (v Vector) Add(a Vector) Vector {
	var (
		data  = v.RawVector().Data
		dataA = a.RawVector().Data
	)
	for i := range data {
		data[i] += dataA[i]
	}
	return v
}

func (v Vector) Subset(i1, i2 int) Vector {
	var (
		data  = v.RawVector().Data
		r     *mat.VecDense
		dataR []float64
	)
	if i1 == i2 {
		i1 = lim(i1, v.Len())
		dataR = make([]float64, 1)
		dataR[0] = data[i1]
	} else {
		i1, i2 = limLoop(i1, i2, v.Len())
		dataR = make([]float64, i2-i1)
		var ind int
		for i := i1; i < i2; i++ {
			dataR[ind] = data[i]
			ind++
		}
	}
	r = mat.NewVecDense(len(dataR), dataR)
	return Vector{V: r, DataP: r.RawVector().Data}
}

func (v Vector) SubsetIndex(I Index) Vector {
	var (
		data  = v.RawVector().Data
		dataR = make([]float64, len(I))
		r     *mat.VecDense
	)
	for i, ind := range I {
		dataR[i] = data[ind]
	}
	r = mat.NewVecDense(len(dataR), dataR)
	return Vector{V: r, DataP: r.RawVector().Data}
}

func (v Vector) Scale(a float64) Vector {
	var (
		data = v.RawVector().Data
	)
	for i := range data {
		data[i] *= a
	}
	return v
}

func (v Vector) AddScalar(a float64) Vector {
	var (
		data = v.RawVector().Data
	)
	for i := range data {
		data[i] += a
	}
	return v
}

func (v Vector) Set(a float64) Vector {
	var (
		data = v.RawVector().Data
	)
	for i := range data {
		data[i] = a
	}
	return v
}

func (v Vector) Apply(f func(float64) float64) Vector {
	var (
		data = v.RawVector().Data
	)
	for i, val := range data {
		data[i] = f(val)
	}
	return v
}

func (v Vector) POW(p int) Vector {
	var (
		data = v.RawVector().Data
	)
	for i, val := range data {
		data[i] = POW(val, p)
	}
	return v
}

func (v Vector) Copy() Vector {
	vv := mat.VecDenseCopyOf(v.V)
	return Vector{
		V:     vv,
		DataP: vv.RawVector().Data,
	}
}

func (v Vector) Find(op EvalOp, target float64, abs bool) (r Vector) {
	var (
		vD = v.RawVector().Data
		rD []float64
	)
	switch op {
	case Equal:
		for i, val := range vD {
			if abs {
				val = math.Abs(val)
			}
			if val == target {
				rD = append(rD, float64(i))
			}
		}
	case Less:
		for i, val := range vD {
			if abs {
				val = math.Abs(val)
			}
			if val < target {
				rD = append(rD, float64(i))
			}
		}
	case Greater:
		for i, val := range vD {
			if abs {
				val = math.Abs(val)
			}
			if val > target {
				rD = append(rD, float64(i))
			}
		}
	case LessOrEqual:
		for i, val := range vD {
			if abs {
				val = math.Abs(val)
			}
			if val <= target {
				rD = append(rD, float64(i))
			}
		}
	case GreaterOrEqual:
		for i, val := range vD {
			if abs {
				val = math.Abs(val)
			}
			if val >= target {
				rD = append(rD, float64(i))
			}
		}
	}
	if len(rD) != 0 {
		r = NewVector(len(rD), rD)
	}
	return
}

func (v Vector) FindOr(op EvalOp, target float64, abs bool, A Vector) (r Vector) {
	var (
		vD = v.RawVector().Data
		vA = A.RawVector().Data
		rD []float64
	)
	switch op {
	case Equal:
		for i, val := range vD {
			t2 := vA[i]
			if abs {
				val = math.Abs(val)
				t2 = math.Abs(t2)
			}
			if val == target || t2 == target {
				rD = append(rD, float64(i))
			}
		}
	case Less:
		for i, val := range vD {
			t2 := vA[i]
			if abs {
				val = math.Abs(val)
				t2 = math.Abs(t2)
			}
			if val < target || t2 < target {
				rD = append(rD, float64(i))
			}
		}
	case Greater:
		for i, val := range vD {
			t2 := vA[i]
			if abs {
				val = math.Abs(val)
				t2 = math.Abs(t2)
			}
			if val > target || t2 > target {
				rD = append(rD, float64(i))
			}
		}
	case LessOrEqual:
		for i, val := range vD {
			t2 := vA[i]
			if abs {
				val = math.Abs(val)
				t2 = math.Abs(t2)
			}
			if val <= target || t2 <= target {
				rD = append(rD, float64(i))
			}
		}
	case GreaterOrEqual:
		for i, val := range vD {
			t2 := vA[i]
			if abs {
				val = math.Abs(val)
				t2 = math.Abs(t2)
			}
			if val >= target || t2 >= target {
				rD = append(rD, float64(i))
			}
		}
	}
	if rD != nil {
		r = NewVector(len(rD), rD)
	}
	return
}

func (v Vector) Concat(w Vector) (r Vector) {
	var (
		v1D = v.RawVector().Data
		v2D = w.RawVector().Data
		N   = len(v1D) + len(v2D)
		rD  = make([]float64, N)
	)
	for i, val := range v1D {
		rD[i] = val
	}
	offset := len(v1D)
	for i, val := range v2D {
		rD[i+offset] = val
	}
	r = NewVector(N, rD)
	return
}

// Non Chainable methods
func (v Vector) Transpose() Matrix {
	var (
		nr, nc = v.V.Dims()
	)
	return NewMatrix(nc, nr, v.RawVector().Data)
}

func (v Vector) ToIndex() (I Index) {
	var (
		data = v.RawVector().Data
	)
	I = make(Index, v.Len())
	for i, val := range data {
		I[i] = int(val)
	}
	return
}

func (v Vector) Min() (min float64) {
	var (
		data = v.RawVector().Data
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
		data = v.RawVector().Data
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
	return NewMatrix(v.V.Len(), 1, v.RawVector().Data)
}

func (v Vector) Mul(w Vector) (A Matrix) {
	var (
		dataV = v.RawVector().Data
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

func (v Vector) Outer(w Vector) (R Matrix) {
	var (
		nr, nc = v.Len(), w.Len()
		dataV  = v.RawVector().Data
		dataW  = w.RawVector().Data
	)
	R = NewMatrix(nr, nc)
	for j := 0; j < nc; j++ {
		x := dataW[j]
		if x == 1 {
			R.SetCol(j, dataV)
		} else {
			Rdata := R.RawMatrix().Data
			// ind = nc * i + j
			for i, val := range dataV {
				Rdata[j+nc*i] = x * val
			}
		}
	}
	return
}

func (v Vector) Minmod() (val float64) {
	return math.Max(0, v.Min())
}

func (v Vector) Print(msgI ...string) (o string) {
	var (
		name = ""
	)
	if len(msgI) != 0 {
		name = msgI[0]
	}
	o = fmt.Sprintf("%s = \n%10.8f\n", name, mat.Formatted(v.V, mat.Squeeze()))
	fmt.Printf(o)
	return
}

func (v Vector) Dot(a Vector) (res float64) {
	var (
		vD, aD = v.Data(), a.Data()
	)
	if len(vD) != len(aD) {
		panic(fmt.Errorf("length of vectors must be the same, have len(v)=%d, len(a)=%d\n", len(vD), len(aD)))
	}
	for i, val := range vD {
		res += val * aD[i]
	}
	return
}

package utils

import (
	"math"

	"gonum.org/v1/gonum/mat"
)

func NewVecRange(min, max int) (V *mat.VecDense) {
	var (
		N int = max - min + 1
		x     = make([]float64, N)
	)
	for i := min; i <= max; i++ {
		x[i] = float64(i)
	}
	V = mat.NewVecDense(N, x)
	return
}

func NewVecConst(N int, val float64) (V *mat.VecDense) {
	var (
		x = make([]float64, N)
	)
	for i := 0; i < N; i++ {
		x[i] = val
	}
	V = mat.NewVecDense(N, x)
	return
}

func VecScalarMult(v mat.Vector, a float64) (vo *mat.VecDense) {
	var (
		d = make([]float64, v.Len())
		N = v.Len()
	)
	for i := 0; i < N; i++ {
		val := v.AtVec(i)
		d[i] = val * a
	}
	return mat.NewVecDense(N, d)
}

func VecScalarAdd(v *mat.VecDense, a float64) (vo *mat.VecDense) {
	var (
		dO = v.RawVector().Data
		d  = make([]float64, v.Len())
		N  = v.Len()
	)
	for i := 0; i < N; i++ {
		d[i] = dO[i] + a
	}
	return mat.NewVecDense(N, d)
}

func VecAbs(v mat.Vector) (vo *mat.VecDense) {
	var (
		d = make([]float64, v.Len())
		N = v.Len()
	)
	for i := 0; i < N; i++ {
		val := v.AtVec(i)
		d[i] = math.Abs(val)
	}
	return mat.NewVecDense(N, d)
}

func VecSquare(v mat.Vector) (vo *mat.VecDense) {
	var (
		d = make([]float64, v.Len())
		N = v.Len()
	)
	for i := 0; i < N; i++ {
		val := v.AtVec(i)
		d[i] = val * val
	}
	return mat.NewVecDense(N, d)
}

func VecSub(V mat.Vector, IndexI interface{}) (R *mat.VecDense) {
	switch Index := IndexI.(type) {
	case mat.Vector:
		R = VecSubV(V, Index)
	case Index:
		var r = make([]float64, len(Index))
		for i, ind := range Index {
			r[i] = V.AtVec(ind)
		}
		R = mat.NewVecDense(len(Index), r)
	}
	return
}

func VecSubV(V, VI mat.Vector) (R *mat.VecDense) {
	// vI should contain a list of indices into v
	var (
		n  = VI.Len()
		r  = make([]float64, n)
		nn = V.Len()
	)
	for i := 0; i < VI.Len(); i++ {
		ival := int(VI.AtVec(i))
		if ival > (nn-1) || ival < 0 {
			return nil
		}
		r[i] = V.AtVec(ival)
	}
	R = mat.NewVecDense(n, r)
	return
}

func VecFind(v *mat.VecDense, op EvalOp, target float64, abs bool) (r *mat.VecDense) {
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
	r = mat.NewVecDense(len(rD), rD)
	return
}

func VecConcat(v1, v2 *mat.VecDense) (r *mat.VecDense) {
	var (
		v1D = v1.RawVector().Data
		v2D = v2.RawVector().Data
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
	r = mat.NewVecDense(N, rD)
	return
}

func VecGetF64(v mat.Vector) (r []float64) {
	r = make([]float64, v.Len())
	for i := 0; i < v.Len(); i++ {
		r[i] = v.AtVec(i)
	}
	return
}

func VecOuter(A, B *mat.VecDense) (R *mat.Dense) {
	var (
		nr, nc = A.Len(), B.Len()
	)
	R = mat.NewDense(nr, nc, nil)
	for j := 0; j < nc; j++ {
		x := B.AtVec(j)
		if x == 1 {
			R.SetCol(j, A.RawVector().Data)
		} else {
			Rdata := R.RawMatrix().Data
			Adata := A.RawVector().Data
			// ind = nc * i + j
			for i, val := range Adata {
				Rdata[j+nc*i] = x * val
			}
		}
	}
	return
}

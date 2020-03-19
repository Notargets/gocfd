package utils

import (
    "gonum.org/v1/gonum/mat"
    "math"
)

func VecRange(min, max int) (V *mat.VecDense){
    var (
        N int = max - min + 1
        x = make([]float64, N)
    )
    for i:=min; i<=max; i++ {
        x[i] = float64(i)
    }
    V = mat.NewVecDense(N, x)
    return
}

func VecConst(val float64, N int) (V *mat.VecDense){
    var (
        x = make([]float64, N)
    )
    for i:=0; i<N; i++ {
        x[i] = val
    }
    V = mat.NewVecDense(N, x)
    return
}

func VecScalarMult(a float64, v mat.Vector) (vo *mat.VecDense) {
    var (
        d = make([]float64, v.Len())
        N = v.Len()
    )
    for i:=0; i<N; i++ {
        val := v.AtVec(i)
        d[i] = val*a
    }
    return mat.NewVecDense(N,d)
}

func VecScalarAdd(a float64, v mat.Vector) (vo *mat.VecDense) {
    var (
        d = make([]float64, v.Len())
        N = v.Len()
    )
    for i:=0; i<N; i++ {
        val := v.AtVec(i)
        d[i] = val+a
    }
    return mat.NewVecDense(N,d)
}

func VecAbs(v mat.Vector) (vo *mat.VecDense) {
    var (
        d = make([]float64, v.Len())
        N = v.Len()
    )
    for i:=0; i<N; i++ {
        val := v.AtVec(i)
        d[i] = math.Abs(val)
    }
    return mat.NewVecDense(N,d)
}

func SquareVector(v mat.Vector) (vo *mat.VecDense) {
    var (
        d = make([]float64, v.Len())
        N = v.Len()
    )
    for i:=0; i<N; i++ {
        val := v.AtVec(i)
        d[i] = val*val
    }
    return mat.NewVecDense(N,d)
}

func SubVector(V, VI mat.Vector) (R *mat.VecDense) {
    // vI should contain a list of indices into v
    var (
        n = VI.Len()
        r = make([]float64, n)
        nn = V.Len()
    )
    for i:=0; i<VI.Len(); i++ {
        ival := int(VI.AtVec(i))
        if ival > (nn-1) || ival < 0 {
            return nil
        }
        r[i] = V.AtVec(ival)
    }
    R = mat.NewVecDense(n, r)
    return
}

func GetFloat64(v *mat.VecDense) (f []float64) {
    f = make([]float64, v.Len())
    for i:=0; i<v.Len(); i++ {
        f[i] = v.AtVec(i)
    }
    return
}


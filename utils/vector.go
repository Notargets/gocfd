package utils

import "gonum.org/v1/gonum/mat"

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

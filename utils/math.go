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

func VecConst(val float64, N int) (V mat.Vector){
    var (
        x = make([]float64, N)
    )
    for i:=0; i<N; i++ {
        x[i] = val
    }
    V = mat.NewVecDense(N, x)
    return
}

func VecScalarMult(a float64, v mat.Vector) (vo mat.Vector) {
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

func NewTriDiagonal(d0, d1, dm1 []float64) (Tri *mat.Dense) {
    dd := make([]float64, len(d0)*len(d0))
    var p1, p2, p3 int
    for j:=0; j<len(d0); j++ {
        for i := 0; i < len(d0); i++ {
            if i == j {
                dd[i+j*len(d0)] = d0[p1]
                p1++
                if i != len(d0)-1 {
                    dd[+1+i+j*len(d0)] = d1[p2]
                    p2++
                }
                if i != 0 {
                    dd[-1 + i+j*len(d0)] = dm1[p3]
                    p3++
                }
            }
        }
    }
    Tri = mat.NewDense(len(d0), len(d0), dd)
    return
}

func NewSymTriDiagonal(d0, d1 []float64) (Tri *mat.SymDense) {
    dd := make([]float64, len(d0)*len(d0))
    var p1, p2 int
    for j:=0; j<len(d0); j++ {
        for i := 0; i < len(d0); i++ {
            if i == j {
                dd[i+j*len(d0)] = d0[p1]
                p1++
                if i != len(d0)-1 {
                    dd[+1+i+j*len(d0)] = d1[p2]
                    p2++
                }
            }
        }
    }
    Tri = mat.NewSymDense(len(d0), dd)
    return
}

func SquareVector(v mat.Vector) (vo mat.Vector) {
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

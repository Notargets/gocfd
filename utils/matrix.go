package utils

import "gonum.org/v1/gonum/mat"

func MatElementInvert(M mat.Matrix) (O *mat.Dense) {
    var (
        d, s = M.Dims()
    )
    O = mat.NewDense(d, s, nil)
    O.CloneFrom(M)
    for j:=0; j<s; j++ {
        for i := 0; i < d; i++ {
            O.Set(i, j, 1./O.At(i, j))
        }
    }
    return
}

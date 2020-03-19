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

func MatSubRow(M *mat.Dense, RowIndices *mat.VecDense) (R *mat.Dense) {
    // RowIndices should contain a list of row indices into M
    var (
        nr, nc = M.Dims()
        nrI = RowIndices.Len()
        rI = RowIndices.RawVector().Data
    )
    R = mat.NewDense(nrI, nc, nil)
    for i, val := range rI {
        valI := int(val)
        if valI > nr-1 || valI < 0 {
            panic("unable to subset row from matrix, index out of bounds")
        }
        rowSlice := M.RawRowView(valI)
        R.SetRow(i, rowSlice)
    }
    return
}

func MatSubCol(M *mat.Dense, ColIndices *mat.VecDense) (R *mat.Dense) {
    // ColIndices should contain a list of row indices into M
    var (
        nr, nc = M.Dims()
        ncI = ColIndices.Len()
        cI = ColIndices.RawVector().Data
    )
    R = mat.NewDense(nr, ncI, nil)
    for i, val := range cI {
        valI := int(val)
        if valI > nc-1 || valI < 0 {
            panic("unable to subset row from matrix, index out of bounds")
        }
        //colSlice := M.ColView(valI)
        R.SetCol(i, VecGetF64(M.ColView(valI)))
    }
    return
}

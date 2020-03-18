package main

import (
    "fmt"
    "github.com/notargets/gophys/DG1D"
    "gonum.org/v1/gonum/mat"
)

const (
    NODETOL = 1.e-12
    K = 10
    N = 8
    Nfp = 1
    NFaces = 2
)
func main() {
    VX, EToV := SimpleMesh1D(0, 2, K)
    fmt.Printf("VX = \n%v\n", mat.Formatted(VX.T(), mat.Squeeze()))
    fmt.Printf("EToV = \n%v\n", mat.Formatted(EToV.T(), mat.Squeeze()))

    J, R, W := DG1D.JacobiGL(0, 0, N)
    fmt.Printf("J = \n%v\n", mat.Formatted(J, mat.Squeeze()))
    fmt.Printf("R = \n%v\n", mat.Formatted(R, mat.Squeeze()))
    fmt.Printf("W = \n%v\n", mat.Formatted(W, mat.Squeeze()))
    V := DG1D.Vandermonde1D(N, R)
    fmt.Printf("V = \n%v\n", mat.Formatted(V, mat.Squeeze()))
    Vinv := mat.NewDense(N+1, N+1, nil)
    if err := Vinv.Inverse(V); err != nil {
        panic("error inverting V")
    }
    fmt.Printf("Vinv = \n%v\n", mat.Formatted(Vinv, mat.Squeeze()))
    Vr := DG1D.GradVandermonde1D(R, N)
    fmt.Printf("Vr = \n%v\n", mat.Formatted(Vr, mat.Squeeze()))
    Dr := mat.NewDense(N+1, N+1, nil)
    Dr.Mul(Vr, Vinv)
    fmt.Printf("Dr = \n%v\n", mat.Formatted(Dr, mat.Squeeze()))
}

func SimpleMesh1D(xmin, xmax float64, K int) (VX *mat.VecDense, EToV *mat.Dense) {
    // K is the number of elements, there are K+1 vertices
    var (
        x = make([]float64, K+1)
        elementVertex = make([]float64, K*2)
    )
    for i:=0; i<K+1; i++ {
        x[i] = (xmax - xmin) * float64(i) / float64(K) + xmin
    }
    var iter int
    for i:=0; i<K; i++ {
        elementVertex[iter] = float64(i)
        elementVertex[iter+1] = float64(i+1)
        iter+=2
    }
    return mat.NewVecDense(K+1, x), mat.NewDense(K, 2, elementVertex)
}



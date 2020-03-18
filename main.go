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
    Np = N+1
)
func main() {
    Startup1D()
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

func Startup1D() {
    VX, EToV := SimpleMesh1D(0, 2, K)

    J, R, W := DG1D.JacobiGL(0, 0, N)
    V := DG1D.Vandermonde1D(N, R)
    Vinv := mat.NewDense(Np, Np, nil)
    if err := Vinv.Inverse(V); err != nil {
        panic("error inverting V")
    }
    Vr := DG1D.GradVandermonde1D(R, N)
    Dr := mat.NewDense(Np, Np, nil)
    Dr.Product(Vr, Vinv)
    LIFT := DG1D.Lift1D(V, Np, NFaces, Nfp)
    NX := DG1D.Normals1D(NFaces, Nfp, K)
    fmt.Printf("LIFT = \n%v\n", mat.Formatted(LIFT, mat.Squeeze()))

    _, _, _, _, _, _ = VX, EToV, J, W, LIFT, NX
}

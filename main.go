package main

import (
    "fmt"
    "gonum.org/v1/gonum/mat"
    "math"
)

func main() {
    var (
        K int
    )
    K = 10
    VX, EToV := SimpleMesh1D(0, 2, K)
    fmt.Printf("VX = \n%v\n", mat.Formatted(VX.T(), mat.Squeeze()))
    fmt.Printf("EToV = \n%v\n", mat.Formatted(EToV.T(), mat.Squeeze()))
    J, x, W, Vr := JacobiGQ(1, 1, 6, false)
    fmt.Printf("J = \n%v\n", mat.Formatted(J, mat.Squeeze()))
    fmt.Printf("x = \n%v\n", mat.Formatted(x, mat.Squeeze()))
    fmt.Printf("Vr = \n%v\n", mat.Formatted(Vr, mat.Squeeze()))
    fmt.Printf("W = \n%v\n", mat.Formatted(W, mat.Squeeze()))
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

/*
func JacobiGL(alpha, beta float64, N int) (X *mat.VecDense){
    var (
        x = make([]float64, N+1)
    )
}
 */

func JacobiGQ(alpha, beta float64, N int, doSort bool) (J *mat.SymDense, X, W mat.Vector, Vr *mat.Dense) {
    var (
        x, w []float64
        ab1, a1, b1, fac float64
        h1, d0, d1 []float64
    )
    if N == 0 {
        x = []float64{-(alpha-beta)/(alpha+beta+2.)}
        w = []float64{2.}
        return nil, mat.NewVecDense(1, x), mat.NewVecDense(1, w), nil
    }

    h1 = make([]float64, N+1)
    for i:=0; i<N+1; i++ {
        h1[i] = 2*float64(i)+alpha+beta
    }

    // main diagonal: diag(-1/2*(alpha^2-beta^2)./(h1+2)./h1)
    d0 = make([]float64, N+1)
    fac = -.5 * (alpha*alpha - beta*beta)
    for i:=0; i<N+1; i++ {
        val := h1[i]
        d0[i] = fac / (val*(val+2.))
    }
    // Handle division by zero
    eps := 1.e-16
    if alpha+beta < 10*eps {
        d0[0] = 0.
    }

    // 1st upper diagonal: diag(2./(h1(1:N)+2).*sqrt((1:N).*((1:N)+alpha+beta) .* ((1:N)+alpha).*((1:N)+beta)./(h1(1:N)+1)./(h1(1:N)+3)),1);
    // for (i=1; i<=N; ++i) { d1(i)=2.0/(h1(i)+2.0)*sqrt(i*(i+alpha+beta)*(i+alpha)*(i+beta)/(h1(i)+1)/(h1(i)+3.0)); }
    var ip1 float64
    d1 = make([]float64, N)
    for i:=0; i<N; i++ {
        ip1 = float64(i+1)
        val := h1[i]
        d1[i] = 2. / (val+2.)
        d1[i] *= math.Sqrt(ip1*(ip1+alpha+beta)*(ip1+alpha)*(ip1+beta)/((val+1.)*(val+3.)))
    }

    J = NewSymTriDiagonal(d0, d1)

    var eig mat.EigenSym
    ok := eig.Factorize(J, true)
    if !ok {
        panic("eigenvalue decomposition failed")
    }
    x = eig.Values(x)
    X = mat.NewVecDense(N+1, x)

    Vr = mat.NewDense(len(x), len(x), nil)
    eig.VectorsTo(Vr)
    W = SquareVector(Vr.RowView(0))
    ab1 = alpha+beta + 1.
    a1 = alpha + 1.
    b1 = beta + 1.
    fac = math.Gamma(a1)*math.Gamma(b1)*math.Pow(2, ab1) / ab1 / math.Gamma(ab1)
    W = VecScalarMult(fac, W)

    return J, X, W, Vr
}

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

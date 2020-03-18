package main

import (
    "fmt"
    "github.com/notargets/gophys/utils"
    "gonum.org/v1/gonum/mat"
    "math"
)

const (
    NODETOL = 1.e-12
)
func main() {
    var (
        K, N, Nfp, NFaces int
    )
    K = 10
    N = 8
    Nfp = 1
    NFaces = 2
    _, _ = Nfp, NFaces
    VX, EToV := SimpleMesh1D(0, 2, K)
    fmt.Printf("VX = \n%v\n", mat.Formatted(VX.T(), mat.Squeeze()))
    fmt.Printf("EToV = \n%v\n", mat.Formatted(EToV.T(), mat.Squeeze()))

    J, Vr, R, W := JacobiGL(0, 0, N)
    fmt.Printf("J = \n%v\n", mat.Formatted(J, mat.Squeeze()))
    fmt.Printf("R = \n%v\n", mat.Formatted(R, mat.Squeeze()))
    fmt.Printf("Vr = \n%v\n", mat.Formatted(Vr, mat.Squeeze()))
    fmt.Printf("W = \n%v\n", mat.Formatted(W, mat.Squeeze()))
    V := Vandermonde1D(N, R)
    fmt.Printf("V = \n%v\n", mat.Formatted(V, mat.Squeeze()))
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

func JacobiGL(alpha, beta float64, N int) (J, Vr mat.Matrix, X, W *mat.VecDense){
    var (
        x = make([]float64, N+1)
        xint mat.Vector
    )
    if N == 1 {
        x[0] = -1
        x[1] = 1
        return nil, nil, mat.NewVecDense(N+1, x), nil
    }
    J, Vr, xint, W = JacobiGQ(alpha + 1, beta + 1, N-2)
    x[0] = -1
    x[N] = 1
    var iter int
    for i:=1; i<N; i++ {
        x[i] = xint.AtVec(iter)
        iter++
    }
    X = mat.NewVecDense(len(x), x)
    return
}

func JacobiGQ(alpha, beta float64, N int) (J, Vr mat.Matrix, X, W *mat.VecDense) {
    var (
        x, w []float64
        fac float64
        h1, d0, d1 []float64
    )
    if N == 0 {
        x = []float64{-(alpha-beta)/(alpha+beta+2.)}
        w = []float64{2.}
        return nil, nil, mat.NewVecDense(1, x), mat.NewVecDense(1, w)
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

    JJ := utils.NewSymTriDiagonal(d0, d1)

    var eig mat.EigenSym
    ok := eig.Factorize(JJ, true)
    if !ok {
        panic("eigenvalue decomposition failed")
    }
    x = eig.Values(x)
    X = mat.NewVecDense(N+1, x)

    VVr := mat.NewDense(len(x), len(x), nil)
    eig.VectorsTo(VVr)
    W = utils.SquareVector(VVr.RowView(0))
    W = utils.VecScalarMult(gamma0(alpha, beta), W)

    return JJ, VVr, X, W
}

func Vandermonde1D(N int, R *mat.VecDense) (V *mat.Dense) {
    V = mat.NewDense(R.Len(), N+1, nil)
    for j:=0; j<N+1; j++ {
        //JP := JacobiP(R, 0, 0, j)
        //for i:=0; i<len(JP); i++ {
        //    fmt.Printf("JP[%d] = \n%v\n", i, JP[i])
        //}
        V.SetCol(j, JacobiP(R, 0, 0, j))
    }
    return
}

func JacobiP(r *mat.VecDense, alpha, beta float64, N int) (p []float64){
    var (
        Nc = r.Len()
    )
    rg := 1./math.Sqrt(gamma0(alpha, beta))
    if N == 0 {
       p = constArray(rg, Nc)
       return
    }
    Np1 := N+1
    pl := make([]float64, Np1*Nc)
    var iter int
    for i:=0; i<Nc; i++ {
        pl[i+iter] = rg
    }

    iter += Nc // Increment to next row
    ab := alpha+beta
    rg1 := 1. / math.Sqrt(gamma1(alpha, beta))
    for i:=0; i<Nc; i++ {
        pl[i+iter] = rg1 * ((ab+2.0)*r.AtVec(i)/2.0 + (alpha-beta)/2.0)
    }

    if N == 1 {
        p = pl[iter:iter+Nc]
        return
    }

    a1 := alpha + 1.
    b1 := beta + 1.
    ab1 := ab + 1.
    aold := 2.0*math.Sqrt(a1*b1/(ab+3.0))/(ab+2.0)
    PL := mat.NewDense(Np1, Nc, pl)
    var xrow []float64
    for i:=0; i<N-1; i++ {
        ip1 := float64(i+1)
        ip2 := float64(ip1+1)
        h1 := 2.0*ip1 + ab
        anew := 2.0/(h1+2.0)*math.Sqrt(ip2*(ip1+ab1)*(ip1+a1)*(ip1+b1)/(h1+1.0)/(h1+3.0))
        bnew := - (alpha*alpha - beta*beta)/h1/(h1+2.0)
        x_bnew := utils.VecConst(-bnew, r.Len())
        x_bnew.AddVec(x_bnew, r)
        xi := PL.RawRowView(i)
        xip1 := PL.RawRowView(i+1)
        xrow = make([]float64, len(xi))
        for j := range xi {
            xrow[j] = 1./(anew*(-aold*xi[j] + x_bnew.AtVec(j)*xip1[j]))
        }
        PL.SetRow(i+2, xrow)
        aold = anew
    }
    p = PL.RawRowView(N)
    return
}

func gamma0(alpha, beta float64) float64 {
    ab1 := alpha+beta + 1.
    a1 := alpha + 1.
    b1 := beta + 1.
    return math.Gamma(a1)*math.Gamma(b1)*math.Pow(2, ab1) / ab1 / math.Gamma(ab1)
}
func gamma1(alpha, beta float64) float64 {
    ab := alpha+beta
    a1 := alpha + 1.
    b1 := beta + 1.
    return a1 * b1 * gamma0(alpha, beta) / (ab+3.0)
}

func constArray(val float64, N int) (v []float64) {
    v = make([]float64, N)
    for i := range v {
        v[i] = val
    }
    return
}

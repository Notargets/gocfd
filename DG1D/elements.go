package DG1D

import (
    "github.com/notargets/gophys/utils"
    "gonum.org/v1/gonum/mat"
    "math"
)

func JacobiGL(alpha, beta float64, N int) (J *mat.SymDense, X, W *mat.VecDense){
    var (
        x = make([]float64, N+1)
        xint mat.Vector
    )
    if N == 1 {
        x[0] = -1
        x[1] = 1
        return nil, mat.NewVecDense(N+1, x), nil
    }
    J, xint, W = JacobiGQ(alpha + 1, beta + 1, N-2)
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

func JacobiGQ(alpha, beta float64, N int) (J *mat.SymDense, X, W *mat.VecDense) {
    var (
        x, w []float64
        fac float64
        h1, d0, d1 []float64
        VVr *mat.Dense
    )
    if N == 0 {
        x = []float64{-(alpha-beta)/(alpha+beta+2.)}
        w = []float64{2.}
        return nil, mat.NewVecDense(1, x), mat.NewVecDense(1, w)
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

    VVr = mat.NewDense(len(x), len(x), nil)
    eig.VectorsTo(VVr)
    W = utils.SquareVector(VVr.RowView(0))
    W = utils.VecScalarMult(gamma0(alpha, beta), W)

    return JJ, X, W
}

func Vandermonde1D(N int, R *mat.VecDense) (V *mat.Dense) {
    V = mat.NewDense(R.Len(), N+1, nil)
    for j:=0; j<N+1; j++ {
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
        p = utils.ConstArray(rg, Nc)
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
            xrow[j] = (-aold*xi[j] + x_bnew.AtVec(j)*xip1[j])/anew
        }
        PL.SetRow(i+2, xrow)
        aold = anew
    }
    p = PL.RawRowView(N)
    return
}

func GradJacobiP(r *mat.VecDense, alpha, beta float64, N int) (p []float64){
    if N == 0 {
        p = make([]float64, r.Len())
        return
    }
    p = JacobiP(r, alpha+1, beta+1, N-1)
    fN := float64(N)
    fac := math.Sqrt(fN*(fN+alpha+beta+1))
    for i, val := range p {
        p[i] = val * fac
    }
    return
}

func GradVandermonde1D(r *mat.VecDense, N int) (Vr *mat.Dense) {
    Vr = mat.NewDense(r.Len(), N+1, nil)
    for i:=0; i<N+1; i++ {
        Vr.SetCol(i, GradJacobiP(r, 0, 0, i))
    }
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

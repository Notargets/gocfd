package DG1D

import (
	"fmt"
	"math"

	"github.com/james-bowman/sparse"
	"github.com/notargets/gocfd/utils"
	"gonum.org/v1/gonum/mat"
)

func JacobiGL(alpha, beta float64, N int) (X, W utils.Vector) {
	var (
		x    = make([]float64, N+1)
		xint utils.Vector
	)
	if N == 1 {
		x[0] = -1
		x[1] = 1
		X = utils.NewVector(N+1, x)
		return
	}
	xint, W = JacobiGQ(alpha+1, beta+1, N-2)
	x[0] = -1
	x[N] = 1
	var iter int
	dataXint := xint.V.RawVector().Data
	for i := 1; i < N; i++ {
		//x[i] = xint.AtVec(iter)
		x[i] = dataXint[iter]
		iter++
	}
	X = utils.NewVector(len(x), x)
	return
}

func JacobiGQ(alpha, beta float64, N int) (X, W utils.Vector) {
	var (
		x, w       []float64
		fac        float64
		h1, d0, d1 []float64
		VVr        *mat.Dense
	)
	if N == 0 {
		x = []float64{-(alpha - beta) / (alpha + beta + 2.)}
		w = []float64{2.}
		return utils.NewVector(len(x), x), utils.NewVector(len(w), w)
	}

	h1 = make([]float64, N+1)
	for i := 0; i < N+1; i++ {
		h1[i] = 2*float64(i) + alpha + beta
	}

	// main diagonal: diag(-1/2*(alpha^2-beta^2)./(h1+2)./h1)
	d0 = make([]float64, N+1)
	fac = -.5 * (alpha*alpha - beta*beta)
	for i := 0; i < N+1; i++ {
		val := h1[i]
		d0[i] = fac / (val * (val + 2.))
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
	for i := 0; i < N; i++ {
		ip1 = float64(i + 1)
		val := h1[i]
		d1[i] = 2. / (val + 2.)
		d1[i] *= math.Sqrt(ip1 * (ip1 + alpha + beta) * (ip1 + alpha) * (ip1 + beta) / ((val + 1.) * (val + 3.)))
	}

	JJ := utils.NewSymTriDiagonal(d0, d1)

	var eig mat.EigenSym
	ok := eig.Factorize(JJ, true)
	if !ok {
		panic("eigenvalue decomposition failed")
	}
	x = eig.Values(x)
	X = utils.NewVector(N+1, x)

	VVr = mat.NewDense(len(x), len(x), nil)
	eig.VectorsTo(VVr)
	W = utils.NewVector(len(x), VVr.RawRowView(0)).POW(2).Scale(gamma0(alpha, beta))
	return X, W
}

func Vandermonde1D(N int, R utils.Vector) (V utils.Matrix) {
	V = utils.NewMatrix(R.Len(), N+1)
	for j := 0; j < N+1; j++ {
		V.SetCol(j, JacobiP(R, 0, 0, j))
	}
	return
}

func JacobiP(r utils.Vector, alpha, beta float64, N int) (p []float64) {
	var (
		Nc = r.Len()
	)
	rg := 1. / math.Sqrt(gamma0(alpha, beta))
	if N == 0 {
		p = utils.ConstArray(rg, Nc)
		return
	}
	Np1 := N + 1
	pl := make([]float64, Np1*Nc)
	var iter int
	for i := 0; i < Nc; i++ {
		pl[i+iter] = rg
	}

	iter += Nc // Increment to next row
	ab := alpha + beta
	rg1 := 1. / math.Sqrt(gamma1(alpha, beta))
	for i := 0; i < Nc; i++ {
		pl[i+iter] = rg1 * ((ab+2.0)*r.AtVec(i)/2.0 + (alpha-beta)/2.0)
	}

	if N == 1 {
		p = pl[iter : iter+Nc]
		return
	}

	a1 := alpha + 1.
	b1 := beta + 1.
	ab1 := ab + 1.
	aold := 2.0 * math.Sqrt(a1*b1/(ab+3.0)) / (ab + 2.0)
	PL := mat.NewDense(Np1, Nc, pl)
	var xrow []float64
	for i := 0; i < N-1; i++ {
		ip1 := float64(i + 1)
		ip2 := float64(ip1 + 1)
		h1 := 2.0*ip1 + ab
		anew := 2.0 / (h1 + 2.0) * math.Sqrt(ip2*(ip1+ab1)*(ip1+a1)*(ip1+b1)/(h1+1.0)/(h1+3.0))
		bnew := -(alpha*alpha - beta*beta) / h1 / (h1 + 2.0)
		x_bnew := utils.NewVector(r.Len()).Set(-bnew)
		x_bnew.Add(r)
		xi := PL.RawRowView(i)
		xip1 := PL.RawRowView(i + 1)
		xrow = make([]float64, len(xi))
		for j := range xi {
			xrow[j] = (-aold*xi[j] + x_bnew.AtVec(j)*xip1[j]) / anew
		}
		PL.SetRow(i+2, xrow)
		aold = anew
	}
	p = PL.RawRowView(N)
	return
}

func GradJacobiP(r utils.Vector, alpha, beta float64, N int) (p []float64) {
	if N == 0 {
		p = make([]float64, r.Len())
		return
	}
	p = JacobiP(r, alpha+1, beta+1, N-1)
	fN := float64(N)
	fac := math.Sqrt(fN * (fN + alpha + beta + 1))
	for i, val := range p {
		p[i] = val * fac
	}
	return
}

func GradVandermonde1D(r utils.Vector, N int) (Vr utils.Matrix) {
	Vr = utils.NewMatrix(r.Len(), N+1)
	for i := 0; i < N+1; i++ {
		Vr.SetCol(i, GradJacobiP(r, 0, 0, i))
	}
	return
}

func Lift1D(V utils.Matrix, Np, Nfaces, Nfp int) (LIFT utils.Matrix) {
	Emat := utils.NewMatrix(Np, Nfaces*Nfp)
	Emat.Set(0, 0, 1)
	Emat.Set(Np-1, 1, 1)
	LIFT = V.Mul(V.Transpose()).Mul(Emat)
	return
}

func Normals1D(Nfaces, Nfp, K int) (NX utils.Matrix) {
	nx := make([]float64, Nfaces*Nfp*K)
	for i := 0; i < K; i++ {
		nx[i] = -1
		nx[i+K] = 1
	}
	NX = utils.NewMatrix(Nfp*Nfaces, K, nx)
	return
}

func GeometricFactors1D(Dr, X utils.Matrix) (J, Rx utils.Matrix) {
	J = Dr.Mul(X)
	Rx = J.Copy().POW(-1)
	return
}

func Connect1D(EToV utils.Matrix) (EToE, EToF utils.Matrix) {
	var (
		NFaces     = 2
		K, _       = EToV.Dims()
		Nv         = K + 1
		TotalFaces = NFaces * K
		vn         = utils.NewVector(2, []float64{0, 1}) // local face to vertex connections
	)
	SpFToV_Tmp := sparse.NewDOK(TotalFaces, Nv)
	var sk int
	for k := 0; k < K; k++ {
		for face := 0; face < NFaces; face++ {
			col := int(vn.AtVec(face))
			SpFToV_Tmp.Set(sk, int(EToV.At(k, col)), 1)
			sk++
		}
	}
	SpFToF := sparse.NewCSR(TotalFaces, TotalFaces, nil, nil, nil)
	SpFToV := SpFToV_Tmp.ToCSR()
	SpFToF.Mul(SpFToV, SpFToV.T())
	for i := 0; i < TotalFaces; i++ {
		v := SpFToF.At(i, i)
		SpFToF.Set(i, i, v-2)
	}
	//fmt.Printf("SpFToV = \n%v\n", mat.Formatted(SpFToV.T(), mat.Squeeze()))
	//fmt.Printf("SpFToF = \n%v\n", mat.Formatted(SpFToF.T(), mat.Squeeze()))
	FacesIndex := utils.MatFind(SpFToF, utils.Equal, 1)
	//faces1, faces2 := utils.MatFind(SpFToF, 1)
	/*
		IVec element1 = floor( (faces1-1)/ Nfaces ) + 1;
		IVec face1    =   mod( (faces1-1), Nfaces ) + 1;
	*/
	element1 := FacesIndex.RI.Apply(func(val int) int { return val / NFaces })
	face1 := FacesIndex.RI.Apply(func(val int) int { return int(math.Mod(float64(val), float64(NFaces))) })
	/*
		IVec element2 = floor( (faces2-1)/ Nfaces ) + 1;
		IVec face2    =   mod( (faces2-1), Nfaces ) + 1;
	*/
	element2 := FacesIndex.CI.Apply(func(val int) int { return val / NFaces })
	face2 := FacesIndex.CI.Apply(func(val int) int { return int(math.Mod(float64(val), float64(NFaces))) })
	/*
	  // Rearrange into Nelements x Nfaces sized arrays
	  IVec ind = sub2ind(K, Nfaces, element1, face1);

	  EToE = outer(Range(1,K), Ones(Nfaces));
	  EToF = outer(Ones(K), Range(1,Nfaces));

	  EToE(ind) = element2;
	  EToF(ind) = face2;
	*/
	EToE = utils.NewRangeOffset(1, K).Outer(utils.NewOnes(NFaces))
	EToF = utils.NewOnes(K).Outer(utils.NewRangeOffset(1, NFaces))
	var I2D utils.Index2D
	var err error
	nr, nc := EToE.Dims()
	if I2D, err = utils.NewIndex2D(nr, nc, element1, face1); err != nil {
		panic(err)
	}
	if err = EToE.IndexedAssign(I2D, element2); err != nil {
		panic(err)
	}
	if err = EToF.IndexedAssign(I2D, face2); err != nil {
		panic(err)
	}
	fmt.Printf("EToE = \n%v\n", mat.Formatted(EToE, mat.Squeeze()))
	fmt.Printf("EToF = \n%v\n", mat.Formatted(EToF, mat.Squeeze()))
	return
}

func BuildMaps1D(VX, FMask utils.Vector,
	X, EToV, EToE, EToF utils.Matrix,
	K, Np, Nfp, NFaces int,
	NODETOL float64) (vmapM, vmapP, mapB, vmapB, mapI, vmapI, mapO, vmapO utils.Index) {
	var (
		NF = Nfp * NFaces
	)
	// number volume nodes consecutively
	nodeids := utils.NewRangeOffset(1, Np*K)

	// find index of face nodes with respect to volume node ordering
	vmapM = utils.NewIndex(Nfp * NFaces * K)
	idsR := utils.NewFromFloat(FMask.RawVector().Data)
	for k1 := 0; k1 < K; k1++ {
		iL1 := k1 * NF
		iL2 := iL1 + NF
		idsL := utils.NewRangeOffset(iL1+1, iL2) // sequential indices for element k1
		if err := vmapM.IndexedAssign(idsL, nodeids.Subset(idsR)); err != nil {
			panic(err)
		}
		idsR.AddInPlace(Np)
	}

	//var one = utils.NewVecConst(Nfp, 1)
	var one = utils.NewVector(Nfp).Set(1)
	vmapP = utils.NewIndex(Nfp * NFaces * K)
	fmt.Printf("X = \n%v\n", mat.Formatted(X, mat.Squeeze()))
	for k1 := 0; k1 < K; k1++ {
		for f1 := 0; f1 < NFaces; f1++ {
			k2 := int(EToE.At(k1, f1))
			f2 := int(EToF.At(k1, f1))
			v1 := int(EToV.At(k1, f1))
			v2 := int(EToV.At(k1, (f1+1)%NFaces))
			refd := math.Sqrt(utils.POW(VX.AtVec(v1)-VX.AtVec(v2), 2))
			skM := k1 * NF
			skP := k2 * NF
			idsM := utils.NewRangeOffset(1+f1*Nfp+skM, (f1+1)*Nfp+skM)
			idsP := utils.NewRangeOffset(1+f2*Nfp+skP, (f2+1)*Nfp+skP)
			vidM := vmapM.Subset(idsM)
			vidP := vmapM.Subset(idsP)
			x1 := X.SubsetVector(vidM)
			x2 := X.SubsetVector(vidP)
			X1 := x1.Outer(one)
			X2 := x2.Outer(one)
			D := X1.Copy().Subtract(X2.Transpose()).POW(2).Apply(math.Sqrt).Apply(math.Abs)
			idMP := utils.MatFind(D, utils.Less, NODETOL*refd)
			idM := idMP.RI
			idP := idMP.CI
			if err := vmapP.IndexedAssign(idM.Add(f1*Nfp+skM), vidP.Subset(idP)); err != nil {
				panic(err)
			}
		}
	}

	// Create list of boundary nodes
	mapB = vmapP.FindVec(utils.Equal, vmapM)
	vmapB = vmapM.Subset(mapB)
	mapI = utils.NewIndex(1)
	mapO = utils.NewIndex(1).Add(K*NFaces - 1)
	vmapI = utils.NewIndex(1)
	vmapO = utils.NewIndex(1).Add(K*Np - 1)
	return
}

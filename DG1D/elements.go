package DG1D

import (
	"fmt"
	"math"

	"github.com/james-bowman/sparse"
	"github.com/notargets/gocfd/utils"
	"gonum.org/v1/gonum/mat"
)

func JacobiGL(alpha, beta float64, N int) (J *mat.SymDense, X, W *mat.VecDense) {
	var (
		x    = make([]float64, N+1)
		xint mat.Vector
	)
	if N == 1 {
		x[0] = -1
		x[1] = 1
		return nil, mat.NewVecDense(N+1, x), nil
	}
	J, xint, W = JacobiGQ(alpha+1, beta+1, N-2)
	x[0] = -1
	x[N] = 1
	var iter int
	for i := 1; i < N; i++ {
		x[i] = xint.AtVec(iter)
		iter++
	}
	X = mat.NewVecDense(len(x), x)
	return
}

func JacobiGQ(alpha, beta float64, N int) (J *mat.SymDense, X, W *mat.VecDense) {
	var (
		x, w       []float64
		fac        float64
		h1, d0, d1 []float64
		VVr        *mat.Dense
	)
	if N == 0 {
		x = []float64{-(alpha - beta) / (alpha + beta + 2.)}
		w = []float64{2.}
		return nil, mat.NewVecDense(1, x), mat.NewVecDense(1, w)
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
	X = mat.NewVecDense(N+1, x)

	VVr = mat.NewDense(len(x), len(x), nil)
	eig.VectorsTo(VVr)
	W = utils.VecSquare(VVr.RowView(0))
	W = utils.VecScalarMult(gamma0(alpha, beta), W)

	return JJ, X, W
}

func Vandermonde1D(N int, R *mat.VecDense) (V *mat.Dense) {
	V = mat.NewDense(R.Len(), N+1, nil)
	for j := 0; j < N+1; j++ {
		V.SetCol(j, JacobiP(R, 0, 0, j))
	}
	return
}

func JacobiP(r *mat.VecDense, alpha, beta float64, N int) (p []float64) {
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
		x_bnew := utils.VecConst(-bnew, r.Len())
		x_bnew.AddVec(x_bnew, r)
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

func GradJacobiP(r *mat.VecDense, alpha, beta float64, N int) (p []float64) {
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

func GradVandermonde1D(r *mat.VecDense, N int) (Vr *mat.Dense) {
	Vr = mat.NewDense(r.Len(), N+1, nil)
	for i := 0; i < N+1; i++ {
		Vr.SetCol(i, GradJacobiP(r, 0, 0, i))
	}
	return
}

func Lift1D(V *mat.Dense, Np, Nfaces, Nfp int) (LIFT *mat.Dense) {
	Emat := mat.NewDense(Np, Nfaces*Nfp, nil)
	Emat.Set(0, 0, 1)
	Emat.Set(Np-1, 1, 1)
	LIFT = mat.NewDense(Np, Nfaces*Nfp, nil)
	LIFT.Product(V, V.T(), Emat)
	return
}

func Normals1D(Nfaces, Nfp, K int) (NX *mat.Dense) {
	nx := make([]float64, Nfaces*Nfp*K)
	for i := 0; i < K; i++ {
		nx[i] = -1
		nx[i+K] = 1
	}
	NX = mat.NewDense(Nfp*Nfaces, K, nx)
	return
}

func GeometricFactors1D(Dr, X *mat.Dense) (J, Rx *mat.Dense) {
	var (
		xd, xs int = X.Dims()
	)
	J = mat.NewDense(xd, xs, nil)
	J.Product(Dr, X)
	Rx = utils.MatElementInvert(J)
	return
}

func Connect1D(EToV *mat.Dense) (EToE, EToF *mat.Dense) {
	var (
		NFaces     = 2
		K, _       = EToV.Dims()
		Nv         = K + 1
		TotalFaces = NFaces * K
		vn         = mat.NewVecDense(2, []float64{0, 1}) // local face to vertex connections
	)
	_, _, _ = Nv, TotalFaces, vn

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
	faces1, faces2 := utils.MatFind(SpFToF, 1)
	/*
		IVec element1 = floor( (faces1-1)/ Nfaces ) + 1;
		IVec face1    =   mod( (faces1-1), Nfaces ) + 1;
	*/
	element1 := faces1.ApplyFunc(func(val int) int { return val / NFaces })
	face1 := faces1.ApplyFunc(func(val int) int { return int(math.Mod(float64(val), float64(NFaces))) })
	/*
		IVec element2 = floor( (faces2-1)/ Nfaces ) + 1;
		IVec face2    =   mod( (faces2-1), Nfaces ) + 1;
	*/
	element2 := faces2.ApplyFunc(func(val int) int { return val / NFaces })
	face2 := faces2.ApplyFunc(func(val int) int { return int(math.Mod(float64(val), float64(NFaces))) })
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
	var err error
	err = utils.MatIndexedAssign(EToE, element1, face1, element2)
	if err != nil {
		panic(err)
	}
	err = utils.MatIndexedAssign(EToF, element1, face1, face2)
	if err != nil {
		panic(err)
	}
	fmt.Printf("EToE = \n%v\n", mat.Formatted(EToE, mat.Squeeze()))
	fmt.Printf("EToF = \n%v\n", mat.Formatted(EToF, mat.Squeeze()))
	return
}

func BuildMaps1D(VX, FMask *mat.VecDense,
	EToV, EToE, EToF *mat.Dense,
	K, Np, Nfp, NFaces int,
	NODETOL float64) (mapM, mapP, vmapM, vmapP, mapB, vmapB *mat.Dense) {
	/*
	   IVec idsL, idsR, idsM,idsP, vidM,vidP, idM,idP;
	   IMat idMP; DMat X1,X2,D;  DVec x1,x2;
	   int k1=0,f1=0, k2=0,f2=0, skP=0, iL1=0,iL2=0;
	   int v1=0, v2=0;  double refd = 0.0;
	   vmapM.resize(Nfp*Nfaces*K); vmapP.resize(Nfp*Nfaces*K);
	   int NF = Nfp*Nfaces;
	*/
	var (
		k2, f2, skP, iL1, iL2, v1, v2 int
		refd                          float64
		NF                            = Nfp * NFaces
		idsL, idsR                    utils.Index
	)
	/*
	   // number volume nodes consecutively
	   IVec nodeids = Range(1,Np*K);
	*/
	nodeids := utils.NewRangeOffset(1, Np*K)
	_, _, _, _, _, _, _, _, _, _ = k2, f2, skP, v1, v2, refd, NF, idsL, idsR, nodeids
	/*
		// find index of face nodes with respect to volume node ordering
		for (k1=1; k1<=K; ++k1) {
			iL1=(k1-1)*NF; iL2=k1*NF;     // define target range in vmapM
			idsL.range(iL1+1, iL2);       // sequential indices for element k1
			idsR = Fmask + (k1-1)*Np;     // offset Fmask for element k1
			vmapM(idsL) = nodeids(idsR);  // map face nodes in element k1
		}
	*/
	for k1 := 0; k1 < K; k1++ {
		iL1 = k1 * NF
		iL2 = iL1 + NF
		idsL = utils.NewRangeOffset(iL1+1, iL2) // sequential indices for element k1
		fmt.Printf("iL1, iL2, idsL = \n%v, %v, %v\n", iL1, iL2, idsL)
		//idsR = FMask + k1*Np                // offset Fmask for element k1
		//vmapM(idsL) = nodeids(idsR)         // map face nodes in element k1
	}
	/*
	   DVec one(Nfp, 1.0);
	   for (k1=1; k1<=K; ++k1) {
	       for (f1=1; f1<=Nfaces; ++f1) {

	           // find neighbor
	           k2 = EToE(k1,f1); f2 = EToF(k1,f1);

	           // reference length of edge
	           v1 = EToV(k1,f1); v2 = EToV(k1, 1+umMOD(f1,Nfaces));
	           refd = sqrt(SQ(VX(v1)-VX(v2)));

	           skM = (k1-1)*NF;  // offset to element k1
	           skP = (k2-1)*NF;  // offset to element k2

	           idsM.range((f1-1)*Nfp+1+skM, f1*Nfp+skM);
	           idsP.range((f2-1)*Nfp+1+skP, f2*Nfp+skP);

	           // find volume node numbers of left and right nodes
	           vidM = vmapM(idsM); vidP = vmapM(idsP);

	           x1 = x(vidM); x2 = x(vidP);
	           X1 = outer(x1,one);
	           X2 = outer(x2,one);

	           // Compute distance matrix
	           D = sqr(X1-trans(X2));

	           idMP = find2D( sqrt(abs(D)), '<', NODETOL*refd);
	           idM=idMP(All,1); idP=idMP(All,2);

	           idM += (f1-1)*Nfp + skM;  // offset ids to {f1,k1}
	           vmapP(idM) = vidP(idP);   // set external element ids

	           idP += (f2-1)*Nfp + skP;  // offset ids to {f2,k2}
	       }
	   }

	   // Create list of boundary nodes
	   mapB = find(vmapP, '=', vmapM);  vmapB = vmapM(mapB);

	   // Inflow and outflow boundaries, single element vectors for this case
	   mapI.resize(1); mapO.resize(1);
	   mapI(1) = 1; mapO(1) = K*Nfaces;
	   vmapI.resize(1); vmapO.resize(1);
	   vmapI(1) = 1; vmapO(1) = K*Np;
	*/
	return
}

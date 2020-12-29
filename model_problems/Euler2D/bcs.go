package Euler2D

import (
	"math"

	"github.com/notargets/gocfd/utils"

	"github.com/notargets/gocfd/model_problems/Euler2D/isentropic_vortex"
)

func (c *Euler) WallBC(k, Kmax int, Q_Face [4]utils.Matrix, ishift int, normal [2]float64, normalFlux [][4]float64) {
	var (
		Nedge = c.dfr.FluxElement.Nedge
		qfD   = Get4DP(Q_Face)
	)
	for i := 0; i < Nedge; i++ {
		ie := i + ishift
		p := c.FS.GetFlowFunctionAtIndex(k+ie*Kmax, qfD, StaticPressure)
		for n := 0; n < 4; n++ {
			normalFlux[i][0] = 0
			normalFlux[i][1] = normal[0] * p
			normalFlux[i][2] = normal[1] * p
			normalFlux[i][3] = 0
		}
	}
}

func (c *Euler) IVortexBC(Time float64, k, Kmax, ishift int, Q_Face [4]utils.Matrix, normal [2]float64) {
	var (
		Nedge   = c.dfr.FluxElement.Nedge
		Nint    = c.dfr.FluxElement.Nint
		qfD     = Get4DP(Q_Face)
		riemann = true
	)
	// Set the flow variables to the exact solution
	X, Y := c.dfr.FluxX.Data(), c.dfr.FluxY.Data()
	for i := 0; i < Nedge; i++ {
		iL := i + ishift
		indFlux := k + (2*Nint+iL)*Kmax
		x, y := X[indFlux], Y[indFlux]
		iv := c.AnalyticSolution.(*isentropic_vortex.IVortex)
		rho, rhoU, rhoV, E := iv.GetStateC(Time, x, y)
		ind := k + iL*Kmax
		var QBC [4]float64
		if riemann {
			QBC = c.RiemannBC(k, Kmax, iL, qfD, [4]float64{rho, rhoU, rhoV, E}, normal)
		} else {
			QBC = [4]float64{rho, rhoU, rhoV, E}
		}
		qfD[0][ind] = QBC[0]
		qfD[1][ind] = QBC[1]
		qfD[2][ind] = QBC[2]
		qfD[3][ind] = QBC[3]
	}
}

func (c *Euler) FarBC(k, Kmax, ishift int, Q_Face [4]utils.Matrix, normal [2]float64) {
	var (
		Nedge = c.dfr.FluxElement.Nedge
		qfD   = Get4DP(Q_Face)
	)
	for i := 0; i < Nedge; i++ {
		iL := i + ishift
		ind := k + iL*Kmax
		QBC := c.RiemannBC(k, Kmax, iL, qfD, c.FS.Qinf, normal)
		qfD[0][ind] = QBC[0]
		qfD[1][ind] = QBC[1]
		qfD[2][ind] = QBC[2]
		qfD[3][ind] = QBC[3]
	}
}

func (c *Euler) RiemannBC(k, Kmax, i int, QQ [4][]float64, QInf [4]float64, normal [2]float64) (Q [4]float64) {
	/*
			Use Riemann invariants along characteristic lines to calculate 1D flow properties normal to boundary
		Rinf = VnormInf - 2 * Cinf / (Gamma -1)
		Rint = VnormInt + 2 * Cint / (Gamma -1)
		Vn = 0.5 * (Rint + Rinf)
		C = 0.25 * (Gamma -1) *(Rint - Rinf)
		Then, project entropy and lateral velocity from the interior, calculate the other primitive variables from that
		P/(rho^Gamma) = constant
	*/
	var (
		qq                 = c.GetQQ(k, Kmax, i, QQ)
		rhoInt, uInt, vInt = qq[0], qq[1] / qq[0], qq[2] / qq[0]
		pInt               = c.FS.GetFlowFunction(qq, StaticPressure)
		CInt               = c.FS.GetFlowFunction(qq, SoundSpeed)
		Gamma              = c.FS.Gamma
		GM1                = Gamma - 1.
		OOGM1              = 1. / GM1
		rhoInf, uInf, vInf = QInf[0], QInf[1] / QInf[0], QInf[2] / QInf[0]
		pInf               = c.FS.Pinf
		CInf               = c.FS.Cinf
		Vtang, Beta        float64
		tangent            = [2]float64{-normal[1], normal[0]}
	)
	VnormInt := normal[0]*uInt + normal[1]*vInt
	VnormInf := normal[0]*uInf + normal[1]*vInf
	Rinf := VnormInf - 2.*CInf*OOGM1
	Rint := VnormInt + 2.*CInt*OOGM1
	Vnorm := 0.5 * (Rint + Rinf)
	C := 0.25 * GM1 * (Rint - Rinf)
	//fmt.Printf("normal = %8.5f,%8.5f\n", normal[0], normal[1])
	switch {
	case VnormInt < 0: // Inflow, entropy and tangent velocity from Qinf
		Vtang = tangent[0]*uInf + tangent[1]*vInf
		Beta = pInf / math.Pow(rhoInf, Gamma)
	case VnormInt >= 0: // Outflow, entropy and tangent velocity from Qint
		Vtang = tangent[0]*uInt + tangent[1]*vInt
		Beta = pInt / math.Pow(rhoInt, Gamma)
	}
	u := Vnorm*normal[0] + Vtang*tangent[0]
	v := Vnorm*normal[1] + Vtang*tangent[1]
	//fmt.Printf("uInt,vInt=%8.5f,%8.5f u,v=%8.5f,%8.5f\n", uInt, vInt, u, v)
	rho := math.Pow(C*C/(Gamma*Beta), OOGM1)
	p := Beta * math.Pow(rho, Gamma)
	Q[0] = rho
	Q[1] = rho * u
	Q[2] = rho * v
	Q[3] = p*OOGM1 + 0.5*rho*(u*u+v*v)
	return
}
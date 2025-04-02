package Euler2D

import (
	"fmt"
	"math"
	"strings"

	"github.com/notargets/gocfd/utils"
)

func (it InitType) String() (txt string) {
	txt = InitPrintNames[it]
	return
}

type FluxType uint

const (
	FLUX_Average FluxType = iota
	FLUX_LaxFriedrichs
	FLUX_Roe
	FLUX_RoeER
)

var (
	FluxNames = map[string]FluxType{
		"average": FLUX_Average,
		"lax":     FLUX_LaxFriedrichs,
		"roe":     FLUX_Roe,
		"roe-er":  FLUX_RoeER,
	}
	FluxPrintNames = []string{"Average", "Lax Friedrichs", "Roe", "Roe-ER"}
)

func (ft FluxType) Print() (txt string) {
	txt = FluxPrintNames[ft]
	return
}

func NewFluxType(label string) (ft FluxType) {
	var (
		ok  bool
		err error
	)
	label = strings.ToLower(label)
	if ft, ok = FluxNames[label]; !ok {
		err = fmt.Errorf("unable to use flux named %s", label)
		panic(err)
	}
	return
}

func (c *Euler) CalculateFluxTransformed(k, Kmax, i int, Jdet, Jinv utils.Matrix, Q [4]utils.Matrix) (Fr, Fs [4]float64) {
	/*
		Note that the transformation applied here closely follows Romero and Jameson, it isn't what it appears in that
		what looks like a broken / wrong chain rule application below is considered in the overall method composition
	*/
	var (
		JdetD = Jdet.DataP[k]
		JinvD = Jinv.DataP[4*k : 4*(k+1)]
	)
	ind := k + i*Kmax
	Fx, Fy := c.CalculateFlux(Q, ind)
	for n := 0; n < 4; n++ {
		Fr[n] = JdetD * (JinvD[0]*Fx[n] + JinvD[1]*Fy[n])
		Fs[n] = JdetD * (JinvD[2]*Fx[n] + JinvD[3]*Fy[n])
	}
	return
}

func (c *Euler) CalculateFlux(QQ [4]utils.Matrix, ind int) (Fx, Fy [4]float64) {
	Fx, Fy = c.FluxCalcMock(QQ[0].DataP[ind], QQ[1].DataP[ind], QQ[2].DataP[ind], QQ[3].DataP[ind])
	return
}

func (c *Euler) FluxCalcBase(rho, rhoU, rhoV, E float64) (Fx, Fy [4]float64) {
	var (
		oorho = 1. / rho
		u     = rhoU * oorho
		v     = rhoV * oorho
		p     = c.FSFar.GetFlowFunctionBase(rho, rhoU, rhoV, E, StaticPressure)
	)
	Fx, Fy =
		[4]float64{rhoU, rhoU*u + p, rhoU * v, u * (E + p)},
		[4]float64{rhoV, rhoV * u, rhoV*v + p, v * (E + p)}
	return
}

func (c *Euler) FluxJacobianTransformed(k, Kmax, i int, Jdet, Jinv utils.Matrix, Q [4]utils.Matrix) (Fr, Gs [16]float64) {
	var (
		JdetD          = Jdet.DataP[k]
		JinvD          = Jinv.DataP[4*k : 4*(k+1)]
		ind            = k + i*Kmax
		q0, q1, q2, q3 = Q[0].DataP[ind], Q[1].DataP[ind], Q[2].DataP[ind], Q[3].DataP[ind]
		Fx, Gy         [16]float64
	)
	Fx, Gy = c.FluxJacobianCalc(q0, q1, q2, q3)
	// fmt.Printf("Fx,Gy[%d] = %v,%v\n", i, Fx, Gy)
	for ii := range Fx { // Transform individual element Flux Jacobian
		Fr[ii] = JdetD * (JinvD[0]*Fx[ii] + JinvD[1]*Gy[ii])
		Gs[ii] = JdetD * (JinvD[2]*Fx[ii] + JinvD[3]*Gy[ii])
	}
	return
}

func (c *Euler) FluxJacobianCalc(rho, rhoU, rhoV, E float64) (Fx, Gy [16]float64) {
	var (
		oorho  = 1. / rho
		u, v   = rhoU * oorho, rhoV * oorho
		u2, v2 = u * u, v * v
		Gamma  = c.FSFar.Gamma
		GM1    = Gamma - 1
		e0     = E * Gamma * oorho
		h0     = (u2 + v2) * GM1
		h1     = h0 - e0
		h2a    = e0 - 0.5*h0
		h2     = h2a - GM1*u2
		h3     = h2a - GM1*v2
	)
	Fx = [16]float64{
		0, 1, 0, 0,
		0.5*h0 - u2, u * (3. - Gamma), -v * GM1, GM1,
		-u * v, v, u, 0,
		u * h1, h2, -u * v * GM1, Gamma * u,
	}
	Gy = [16]float64{
		0, 0, 1, 0,
		-u * v, v, u, 0,
		0.5*h0 - v2, -u * GM1, v * (3. - Gamma), GM1,
		v * h1, -u * v * GM1, h3, Gamma * v,
	}
	return
}

func (c *Euler) AvgFlux(kL, kR, KmaxL, KmaxR, shiftL, shiftR int,
	Q_FaceL, Q_FaceR [4]utils.Matrix, normal [2]float64, normalFlux [][4]float64) {
	var (
		Nedge = c.DFR.FluxElement.NpEdge
	)
	averageFluxN := func(fx1, fy1, fx2, fy2 [4]float64, normal [2]float64) (fnorm [4]float64) {
		var (
			fave [2][4]float64
		)
		for n := 0; n < 4; n++ {
			fave[0][n] = 0.5 * (fx1[n] + fx2[n])
			fave[1][n] = 0.5 * (fy1[n] + fy2[n])
			fnorm[n] = normal[0]*fave[0][n] + normal[1]*fave[1][n]
		}
		return
	}
	for i := 0; i < Nedge; i++ {
		iL := i + shiftL
		iR := Nedge - 1 - i + shiftR // Shared edges run in reverse order relative to each other
		indL, indR := kL+iL*KmaxL, kR+iR*KmaxR
		FxL, FyL := c.CalculateFlux(Q_FaceL, indL)
		FxR, FyR := c.CalculateFlux(Q_FaceR, indR) // Reverse the right edge to match
		normalFlux[i] = averageFluxN(FxL, FyL, FxR, FyR, normal)
	}
}

func (c *Euler) LaxFlux(kL, kR, KmaxL, KmaxR, shiftL, shiftR int,
	Q_FaceL, Q_FaceR [4]utils.Matrix, normal [2]float64, normalFlux [][4]float64) {
	var (
		Nedge                                  = c.DFR.FluxElement.NpEdge
		EL, rhoL, rhoUL, rhoVL, uL, vL, pL, CL float64
		ER, rhoR, rhoUR, rhoVR, uR, vR, pR, CR float64
		nx, ny                                 = normal[0], normal[1]
	)
	for i := 0; i < Nedge; i++ {
		iL := i + shiftL
		iR := Nedge - 1 - i + shiftR // Shared edges run in reverse order relative to each other
		indL, indR := kL+iL*KmaxL, kR+iR*KmaxR
		rhoL, rhoR = Q_FaceL[0].DataP[indL], Q_FaceR[0].DataP[indR]
		rhoUL, rhoVL = Q_FaceL[1].DataP[indL], Q_FaceL[2].DataP[indL]
		rhoUR, rhoVR = Q_FaceR[1].DataP[indR], Q_FaceR[2].DataP[indR]
		EL, ER = Q_FaceL[3].DataP[indL], Q_FaceR[3].DataP[indR]
		uL, vL = rhoUL/rhoL, rhoVL/rhoL
		uR, vR = rhoUR/rhoR, rhoVR/rhoR
		pL, pR = c.FSFar.GetFlowFunction(Q_FaceL, indL, StaticPressure), c.FSFar.GetFlowFunction(Q_FaceR, indR, StaticPressure)
		CL, CR = c.FSFar.GetFlowFunction(Q_FaceL, indL, SoundSpeed), c.FSFar.GetFlowFunction(Q_FaceR, indR, SoundSpeed)
		maxV := math.Max(math.Sqrt(uL*uL+vL*vL)+CL, math.Sqrt(uR*uR+vR*vR)+CR)
		normalFlux[i][0] = 0.5 * (nx*(rhoUL+rhoUR) + ny*(rhoVL+rhoVR))
		normalFlux[i][1] = 0.5 * (nx*(rhoUL*uL+rhoUR*uR+pL+pR) + ny*(rhoUL*vL+rhoUR*vR))
		normalFlux[i][2] = 0.5 * (nx*(rhoVL*uL+rhoVR*uR) + ny*(rhoVL*vL+rhoVR*vR+pL+pR))
		normalFlux[i][3] = 0.5 * (nx*((pL+EL)*uL+(pR+ER)*uR) + ny*((pL+EL)*vL+(pR+ER)*vR))
		for n := 0; n < 4; n++ {
			normalFlux[i][n] += 0.5 * maxV * (Q_FaceL[n].DataP[indL] - Q_FaceR[n].DataP[indR])
		}
	}
}

func (c *Euler) LaxFlux2(kL, kR, KmaxL, KmaxR, shiftL, shiftR int,
	Q_FaceL, Q_FaceR [4]utils.Matrix, Flux_FaceL, Flux_FaceR [2][4]utils.Matrix,
	normal [2]float64, normalFlux [][4]float64) {
	var (
		Nedge                          = c.DFR.FluxElement.NpEdge
		rhoL, rhoUL, rhoVL, uL, vL, CL float64
		rhoR, rhoUR, rhoVR, uR, vR, CR float64
		nx, ny                         = normal[0], normal[1]
	)
	for i := 0; i < Nedge; i++ {
		iL := i + shiftL
		iR := Nedge - 1 - i + shiftR // Shared edges run in reverse order relative to each other
		indL, indR := kL+iL*KmaxL, kR+iR*KmaxR
		rhoL, rhoR = Q_FaceL[0].DataP[indL], Q_FaceR[0].DataP[indR]
		rhoUL, rhoVL = Q_FaceL[1].DataP[indL], Q_FaceL[2].DataP[indL]
		rhoUR, rhoVR = Q_FaceR[1].DataP[indR], Q_FaceR[2].DataP[indR]
		uL, vL = rhoUL/rhoL, rhoVL/rhoL
		uR, vR = rhoUR/rhoR, rhoVR/rhoR
		CL, CR = c.FSFar.GetFlowFunction(Q_FaceL, indL, SoundSpeed), c.FSFar.GetFlowFunction(Q_FaceR, indR, SoundSpeed)
		maxV := math.Max(math.Sqrt(uL*uL+vL*vL)+CL, math.Sqrt(uR*uR+vR*vR)+CR)
		for n := 0; n < 4; n++ {
			nL := nx*Flux_FaceL[0][n].DataP[indL] + ny*Flux_FaceL[1][n].DataP[indL]
			nR := nx*Flux_FaceR[0][n].DataP[indR] + ny*Flux_FaceR[1][n].DataP[indR]
			normalFlux[i][n] = 0.5 * (nL + nR + maxV*(Q_FaceL[n].DataP[indL]-Q_FaceR[n].DataP[indR]))
		}
	}
}

func (c *Euler) RoeFlux(kL, kR, KmaxL, KmaxR, shiftL, shiftR int,
	Q_FaceL, Q_FaceR [4]utils.Matrix, normal [2]float64, normalFlux [][4]float64) {
	// fmt.Printf("here 1\n")
	var (
		Nedge            = c.DFR.FluxElement.NpEdge
		rhoL, uL, vL, pL float64
		rhoR, uR, vR, pR float64
		hL, hR           float64
		Gamma            = c.FSFar.Gamma
		GM1              = Gamma - 1
	)
	rotate := func(rhoU, rhoV, nx, ny float64) (rhoUr, rhoVr float64) {
		rhoUr = rhoU*nx + rhoV*ny
		rhoVr = rhoU*(-ny) + rhoV*nx
		return
	}
	for i := 0; i < Nedge; i++ {
		iL := i + shiftL
		iR := Nedge - 1 - i + shiftR // Shared edges run in reverse order relative to each other
		indL, indR := kL+iL*KmaxL, kR+iR*KmaxR
		// Rotate the momentum into face normal coordinates before calculating fluxes
		/*
			rhoULr, rhoVLr := Q_FaceL[1].DataP[indL]*normal[0]+Q_FaceL[2].DataP[indL]*normal[1],
				-Q_FaceL[1].DataP[indL]*normal[1]+Q_FaceL[2].DataP[indL]*normal[0]
			rhoURr, rhoVRr := Q_FaceR[1].DataP[indR]*normal[0]+Q_FaceR[2].DataP[indR]*normal[1],
				-Q_FaceR[1].DataP[indR]*normal[1]+Q_FaceR[2].DataP[indR]*normal[0]
		*/
		rhoULr, rhoVLr := rotate(Q_FaceL[1].DataP[indL], Q_FaceL[2].DataP[indL], normal[0], normal[1])
		rhoURr, rhoVRr := rotate(Q_FaceR[1].DataP[indR], Q_FaceR[2].DataP[indR], normal[0], normal[1])
		rhoL, uL, vL = Q_FaceL[0].DataP[indL], rhoULr/Q_FaceL[0].DataP[indL], rhoVLr/Q_FaceL[0].DataP[indL]
		rhoR, uR, vR = Q_FaceR[0].DataP[indR], rhoURr/Q_FaceR[0].DataP[indR], rhoVRr/Q_FaceR[0].DataP[indR]
		// rhoL, uL, vL = Q_FaceL[0].DataP[indL], Q_FaceL[1].DataP[indL]/Q_FaceL[0].DataP[indL], Q_FaceL[2].DataP[indL]/Q_FaceL[0].DataP[indL]
		// rhoR, uR, vR = Q_FaceR[0].DataP[indR], Q_FaceR[1].DataP[indR]/Q_FaceR[0].DataP[indR], Q_FaceR[2].DataP[indR]/Q_FaceR[0].DataP[indR]
		pL, pR = c.FSFar.GetFlowFunction(Q_FaceL, indL, StaticPressure), c.FSFar.GetFlowFunction(Q_FaceR, indR, StaticPressure)
		/*
		   HM = (EnerM+pM).dd(rhoM);  HP = (EnerP+pP).dd(rhoP);
		*/
		// Enthalpy
		hL, hR = (Q_FaceL[3].DataP[indL]+pL)/rhoL, (Q_FaceR[3].DataP[indR]+pR)/rhoR
		/*
			rhoMs = sqrt(rhoM); rhoPs = sqrt(rhoP);
			rhoMsPs = rhoMs + rhoPs;

			rho = rhoMs.dm(rhoPs);
			u   = (rhoMs.dm(uM) + rhoPs.dm(uP)).dd(rhoMsPs);
			v   = (rhoMs.dm(vM) + rhoPs.dm(vP)).dd(rhoMsPs);
			H   = (rhoMs.dm(HM) + rhoPs.dm(HP)).dd(rhoMsPs);
			c2 = gm1 * (H - 0.5*(sqr(u)+sqr(v)));
			C = sqrt(c2);
		*/
		// Compute Roe average variables
		rhoLs, rhoRs := math.Sqrt(rhoL), math.Sqrt(rhoR)
		rhoLsRs := rhoLs + rhoRs

		rho := rhoLs * rhoRs
		u := (rhoLs*uL + rhoRs*uR) / rhoLsRs
		v := (rhoLs*vL + rhoRs*vR) / rhoLsRs
		h := (rhoLs*hL + rhoRs*hR) / rhoLsRs
		c2 := GM1 * (h - 0.5*(u*u+v*v))
		C := math.Sqrt(c2)
		/*
		   dW1 = -0.5*(rho.dm(uP-uM)).dd(C) + 0.5*(pP-pM).dd(c2);
		   dW2 = (rhoP-rhoM) - (pP-pM).dd(c2);
		   dW3 = rho.dm(vP-vM);
		   dW4 = 0.5*(rho.dm(uP-uM)).dd(C) + 0.5*(pP-pM).dd(c2);

		   dW1 = abs(u-C).dm(dW1);
		   dW2 = abs(u  ).dm(dW2);
		   dW3 = abs(u  ).dm(dW3);
		   dW4 = abs(u+C).dm(dW4);
		*/
		// Riemann fluxes
		dW1 := -0.5*(rho*(uR-uL))/C + 0.5*(pR-pL)/c2
		dW2 := (rhoR - rhoL) - (pR-pL)/c2
		dW3 := rho * (vR - vL)
		dW4 := 0.5*(rho*(uR-uL))/C + 0.5*(pR-pL)/c2
		dW1 = math.Abs(u-C) * dW1
		dW2 = math.Abs(u) * dW2
		dW3 = math.Abs(u) * dW3
		dW4 = math.Abs(u+C) * dW4
		/*
		   DMat fx = (fxQP+fxQM)/2.0;
		   fx(All,1) -= (dW1               + dW2                                   + dW4              )/2.0;
		   fx(All,2) -= (dW1.dm(u-C)       + dW2.dm(u)                             + dW4.dm(u+C)      )/2.0;
		   fx(All,3) -= (dW1.dm(v)         + dW2.dm(v)                 + dW3       + dW4.dm(v)        )/2.0;
		   fx(All,4) -= (dW1.dm(H-u.dm(C)) + dW2.dm(sqr(u)+sqr(v))/2.0 + dW3.dm(v) + dW4.dm(H+u.dm(C)))/2.0;
		*/
		// Form Roe FluxIndex
		// Ave of normal component of flux
		// TODO: This is where the direct values of the Left and Right RT face
		// TODO: fluxes should be used. We should sub the F_L and F_R here
		// TODO: The dissipation term below is computed as given from the
		// TODO: scalars
		normalFlux[i][0] = 0.5 * (rhoULr + rhoURr)
		normalFlux[i][1] = 0.5 * (rhoULr*uL + rhoURr*uR + +pL + pR)
		normalFlux[i][2] = 0.5 * (rhoVLr*uL + rhoVRr*uR)
		normalFlux[i][3] = 0.5 * ((pL+Q_FaceL[3].DataP[indL])*uL + (pR+Q_FaceR[3].DataP[indR])*uR)

		normalFlux[i][0] -= 0.5 * (dW1 + dW2 + dW4)
		normalFlux[i][1] -= 0.5 * (dW1*(u-C) + dW2*u + dW4*(u+C))
		normalFlux[i][2] -= 0.5 * (dW1*v + dW2*v + dW3 + dW4*v)
		normalFlux[i][3] -= 0.5 * (dW1*(h-u*C) + 0.5*dW2*(u*u+v*v) + dW3*v + dW4*(h+u*C))
		/*
		   flux = fx;    fx2.borrow(Ngf, fx.pCol(2)); fx3.borrow(Ngf, fx.pCol(3));
		   flux(All,2) = lnx.dm(fx2) - lny.dm(fx3);
		   flux(All,3) = lny.dm(fx2) + lnx.dm(fx3);
		*/
		// rotate back to Cartesian
		normalFlux[i][1], normalFlux[i][2] = normal[0]*normalFlux[i][1]-normal[1]*normalFlux[i][2],
			normal[1]*normalFlux[i][1]+normal[0]*normalFlux[i][2]
	}
}

func (c *Euler) RoeERFlux(kL, kR, KmaxL, KmaxR, shiftL, shiftR int,
	Q_FaceL, Q_FaceR [4]utils.Matrix, normal [2]float64, normalFlux [][4]float64) {
	var (
		Nedge                          = c.DFR.FluxElement.NpEdge
		rhoL, uL, vL, pL, EL, HL, UL   float64
		rhoR, uR, vR, pR, ER, HR, UR   float64
		dU, dP, dRho, dRhoU, dRhoV, dE float64
		Gamma                          = c.FSFar.Gamma
		GM1                            = Gamma - 1
		nx, ny                         = normal[0], normal[1]
		rho, u, v, h, H, U, C          float64
		sqrt, abs, min, max, sign      = math.Sqrt, math.Abs, math.Min, math.Max, math.Copysign
	)
	/*
		Need:
		- QL and QR using Enthalpy for the:
			- Fbar average
			- Upwind term with delta Q
		- UL, PL and UR, PR for the Fbar average - U is the face normal velocity
		- QRoe for deltaUu term in Fd
		- URoe, CRoe for deltaPu and deltaPp terms and coefficients
	*/
	for i := 0; i < Nedge; i++ {
		iL := i + shiftL
		iR := Nedge - 1 - i + shiftR // Shared edges run in reverse order relative to each other
		indL, indR := kL+iL*KmaxL, kR+iR*KmaxR
		rhoL, rhoR = Q_FaceL[0].DataP[indL], Q_FaceR[0].DataP[indR]
		ooRhoL, ooRhoR := 1./rhoL, 1./rhoR
		rhoLs, rhoRs := sqrt(rhoL), sqrt(rhoR)
		uL, vL = Q_FaceL[1].DataP[indL]*ooRhoL, Q_FaceL[2].DataP[indL]*ooRhoL
		uR, vR = Q_FaceR[1].DataP[indR]*ooRhoR, Q_FaceR[2].DataP[indR]*ooRhoR
		EL, ER = Q_FaceL[3].DataP[indL], Q_FaceR[3].DataP[indR]
		UL, UR = nx*uL+ny*vL, nx*uR+ny*vR // Calculate face normal velocity left and right
		pL, pR = c.FSFar.GetFlowFunction(Q_FaceL, indL, StaticPressure), c.FSFar.GetFlowFunction(Q_FaceR, indR, StaticPressure)
		HL, HR = EL+pL, ER+pR
		// Roe averaged variables
		ooRs := 1. / (rhoLs + rhoRs)
		u, v, h = (rhoLs*uL+rhoRs*uR)*ooRs, (rhoLs*vL+rhoRs*vR)*ooRs, (HL+HR)*ooRs
		rho = rhoLs * rhoRs
		H = h * rho
		U = nx*u + ny*v
		C2 := GM1 * (h - 0.5*(u*u+v*v))
		C = sqrt(C2)
		ooC := 1. / C
		Uabs := abs(U)
		/*
			Calculate Fbar
		*/
		normalFlux[i] = [4]float64{
			0.5 * (UL*rhoL + UR*rhoR),
			0.5 * (UL*rhoL*uL + pL*nx + UR*rhoR*uR + pR*nx),
			0.5 * (UL*rhoL*vL + pL*ny + UR*rhoR*vR + pR*ny),
			0.5 * (UL*HL + UR*HR)}

		/*
			Calculate Fd
		*/
		Uef := 0.05 * C // 0.05 is "epsilon", a constant in Roe-ER
		du, dv := (uR - uL), (vR - vL)
		deltaV2 := du*du + dv*dv
		ooVmag := 1. / sqrt(u*u+v*v)
		var n1x, n1y, n2x, n2y float64
		if deltaV2 < 0.01*C2 { // the number "0.01" is (delta=0.1)^2, the constant "delta" in Roe-ER
			n1x, n1y = nx, ny
		} else {
			n1x, n1y = ooVmag*du, ooVmag*dv
		}
		// Sigma parameter: calculate n2 = (n1 x n) x n1
		n2x, n2y = n1y*(nx*n1y-n1x*ny), -n1x*(nx*n1y-n1x*ny)
		alp1, alp2 := nx*n1x+ny*n1y, nx*n2x+ny*n2y
		// Get rotated velocities
		U1x, U1y := n1x*u, n1y*v
		U2x, U2y := n2x*u, n2y*v
		Urot := sqrt(alp1*alp1*(U1x*U1x+U1y*U1y)) + sqrt(alp2*alp2*(U2x*U2x+U2y*U2y))
		sigma := max(Uabs, min(Uef, Urot))

		UabsPrime := Uabs - 0.25*max(0, UR-UL)*(sign(U+C, 1)-sign(U-C, 1))
		dU, dP, dRho, dRhoU, dRhoV, dE = UR-UL, pR-pL, rhoR-rhoL, rhoR*uR-rhoL*uL, rhoR*vR-rhoL*vL, ER-EL
		dPu := rho * dU * max(0, C-UabsPrime)
		swt := sign(U, 1) * min(UabsPrime, C)
		dPp := swt * dP * ooC
		dUu := swt * dU * ooC

		normalFlux[i][0] -= 0.5 * (sigma*dRho + (dPu+dPp)*0 + dUu*rho)
		normalFlux[i][1] -= 0.5 * (sigma*dRhoU + (dPu+dPp)*nx + dUu*rho*u)
		normalFlux[i][2] -= 0.5 * (sigma*dRhoV + (dPu+dPp)*ny + dUu*rho*v)
		normalFlux[i][3] -= 0.5 * (sigma*dE + (dPu+dPp)*ny + dUu*H)
	}
}

func (c *Euler) RoeFlux2(kL, kR, KmaxL, KmaxR, shiftL, shiftR int,
	Q_FaceL, Q_FaceR [4]utils.Matrix, Flux_FaceL, Flux_FaceR [2][4]utils.Matrix,
	normal [2]float64, normalFlux [][4]float64) {
	// fmt.Printf("here 2\n")
	var (
		Nedge            = c.DFR.FluxElement.NpEdge
		rhoL, uL, vL, pL float64
		rhoR, uR, vR, pR float64
		hL, hR           float64
		Gamma            = c.FSFar.Gamma
		GM1              = Gamma - 1
	)
	rotate := func(rhoU, rhoV, nx, ny float64) (rhoUr, rhoVr float64) {
		rhoUr = rhoU*nx + rhoV*ny
		rhoVr = rhoU*(-ny) + rhoV*nx
		return
	}
	for i := 0; i < Nedge; i++ {
		iL := i + shiftL
		iR := Nedge - 1 - i + shiftR // Shared edges run in reverse order relative to each other
		indL, indR := kL+iL*KmaxL, kR+iR*KmaxR
		// Rotate the momentum into face normal coordinates before calculating fluxes
		rhoULr, rhoVLr := rotate(Q_FaceL[1].DataP[indL], Q_FaceL[2].DataP[indL], normal[0], normal[1])
		rhoURr, rhoVRr := rotate(Q_FaceR[1].DataP[indR], Q_FaceR[2].DataP[indR], normal[0], normal[1])
		rhoL, uL, vL = Q_FaceL[0].DataP[indL], rhoULr/Q_FaceL[0].DataP[indL], rhoVLr/Q_FaceL[0].DataP[indL]
		rhoR, uR, vR = Q_FaceR[0].DataP[indR], rhoURr/Q_FaceR[0].DataP[indR], rhoVRr/Q_FaceR[0].DataP[indR]
		pL, pR = c.FSFar.GetFlowFunction(Q_FaceL, indL, StaticPressure), c.FSFar.GetFlowFunction(Q_FaceR, indR, StaticPressure)
		/*
		   HM = (EnerM+pM).dd(rhoM);  HP = (EnerP+pP).dd(rhoP);
		*/
		// Enthalpy
		hL, hR = (Q_FaceL[3].DataP[indL]+pL)/rhoL, (Q_FaceR[3].DataP[indR]+pR)/rhoR
		// Compute Roe average variables
		rhoLs, rhoRs := math.Sqrt(rhoL), math.Sqrt(rhoR)
		rhoLsRs := rhoLs + rhoRs

		rho := rhoLs * rhoRs
		u := (rhoLs*uL + rhoRs*uR) / rhoLsRs
		v := (rhoLs*vL + rhoRs*vR) / rhoLsRs
		h := (rhoLs*hL + rhoRs*hR) / rhoLsRs
		c2 := GM1 * (h - 0.5*(u*u+v*v))
		C := math.Sqrt(c2)
		// Riemann fluxes
		dW1 := -0.5*(rho*(uR-uL))/C + 0.5*(pR-pL)/c2
		dW2 := (rhoR - rhoL) - (pR-pL)/c2
		dW3 := rho * (vR - vL)
		dW4 := 0.5*(rho*(uR-uL))/C + 0.5*(pR-pL)/c2
		dW1 = math.Abs(u-C) * dW1
		dW2 = math.Abs(u) * dW2
		dW3 = math.Abs(u) * dW3
		dW4 = math.Abs(u+C) * dW4
		normalFlux[i][0] = -0.5 * (dW1 + dW2 + dW4)
		normalFlux[i][1] = -0.5 * (dW1*(u-C) + dW2*u + dW4*(u+C))
		normalFlux[i][2] = -0.5 * (dW1*v + dW2*v + dW3 + dW4*v)
		normalFlux[i][3] = -0.5 * (dW1*(h-u*C) + 0.5*dW2*(u*u+v*v) + dW3*v + dW4*(h+u*C))

		// rotate back to Cartesian
		normalFlux[i][1], normalFlux[i][2] = normal[0]*normalFlux[i][1]-normal[1]*normalFlux[i][2],
			normal[1]*normalFlux[i][1]+normal[0]*normalFlux[i][2]

		for n := 0; n < 4; n++ {
			nL := normal[0]*Flux_FaceL[0][n].DataP[indL] + normal[1]*Flux_FaceL[1][n].DataP[indL]
			nR := normal[0]*Flux_FaceR[0][n].DataP[indR] + normal[1]*Flux_FaceR[1][n].DataP[indR]
			normalFlux[i][n] += 0.5 * (nL + nR)
			// _, _ = nL, nR
		}
	}
}

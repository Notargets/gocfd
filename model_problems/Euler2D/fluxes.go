package Euler2D

import (
	"fmt"
	"math"
	"strings"

	"github.com/notargets/gocfd/utils"
)

func (it InitType) Print() (txt string) {
	txt = InitPrintNames[it]
	return
}

type FluxType uint

const (
	FLUX_Average FluxType = iota
	FLUX_LaxFriedrichs
	FLUX_Roe
)

var (
	FluxNames = map[string]FluxType{
		"average": FLUX_Average,
		"lax":     FLUX_LaxFriedrichs,
		"roe":     FLUX_Roe,
	}
	FluxPrintNames = []string{"Average", "Lax Friedrichs", "Roe"}
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

func (c *Euler) CalculateFluxTransformed(k, Kmax, i int, Jdet, Jinv utils.Matrix, QQ [4][]float64) (Fr, Fs [4]float64) {
	var (
		JdetD = Jdet.Data()[k]
		JinvD = Jinv.Data()[4*k : 4*(k+1)]
	)
	Fx, Fy := c.CalculateFlux(k, Kmax, i, QQ)
	for n := 0; n < 4; n++ {
		Fr[n] = JdetD * (JinvD[0]*Fx[n] + JinvD[1]*Fy[n])
		Fs[n] = JdetD * (JinvD[2]*Fx[n] + JinvD[3]*Fy[n])
	}
	return
}

func (c *Euler) CalculateFlux(k, Kmax, i int, QQ [4][]float64) (Fx, Fy [4]float64) {
	// From https://www.theoretical-physics.net/dev/fluid-dynamics/euler.html
	var (
		ind = k + i*Kmax
		qq  = [4]float64{QQ[0][ind], QQ[1][ind], QQ[2][ind], QQ[3][ind]}
	)
	Fx, Fy = c.FluxCalcMock(qq)
	return
}

func (c *Euler) FluxCalc(q [4]float64) (Fx, Fy [4]float64) {
	var (
		rho, rhoU, rhoV, E = q[0], q[1], q[2], q[3]
		oorho              = 1. / rho
		u                  = rhoU * oorho
		v                  = rhoV * oorho
		p                  = c.FS.GetFlowFunction(q, StaticPressure)
	)
	Fx, Fy =
		[4]float64{rhoU, rhoU*u + p, rhoU * v, u * (E + p)},
		[4]float64{rhoV, rhoV * u, rhoV*v + p, v * (E + p)}
	return
}

func (c *Euler) AvgFlux(kL, kR, KmaxL, KmaxR, shiftL, shiftR int,
	Q_FaceL, Q_FaceR [4]utils.Matrix, normal [2]float64, normalFlux, normalFluxReversed [][4]float64) {
	var (
		Nedge      = c.dfr.FluxElement.Nedge
		qfDL, qfDR = Get4DP(Q_FaceL), Get4DP(Q_FaceR)
	)
	averageFluxN := func(fx1, fy1, fx2, fy2 [4]float64, normal [2]float64) (fnorm [4]float64, fnormR [4]float64) {
		var (
			fave [2][4]float64
		)
		for n := 0; n < 4; n++ {
			fave[0][n] = 0.5 * (fx1[n] + fx2[n])
			fave[1][n] = 0.5 * (fy1[n] + fy2[n])
			fnorm[n] = normal[0]*fave[0][n] + normal[1]*fave[1][n]
			fnormR[n] = -fnorm[n]
		}
		return
	}
	for i := 0; i < Nedge; i++ {
		iL := i + shiftL
		iR := Nedge - 1 - i + shiftR // Shared edges run in reverse order relative to each other
		FxL, FyL := c.CalculateFlux(kL, KmaxL, iL, qfDL)
		FxR, FyR := c.CalculateFlux(kR, KmaxR, iR, qfDR) // Reverse the right edge to match
		normalFlux[i], normalFluxReversed[Nedge-1-i] = averageFluxN(FxL, FyL, FxR, FyR, normal)
	}
}

func (c *Euler) LaxFlux(kL, kR, KmaxL, KmaxR, shiftL, shiftR int,
	Q_FaceL, Q_FaceR [4]utils.Matrix, normal [2]float64, normalFlux, normalFluxReversed [][4]float64) {
	var (
		Nedge                = c.dfr.FluxElement.Nedge
		rhoL, uL, vL, pL, CL float64
		rhoR, uR, vR, pR, CR float64
		qfDL, qfDR           = Get4DP(Q_FaceL), Get4DP(Q_FaceR)
	)
	maxVF := func(u, v, p, rho, C float64) (vmax float64) {
		vmax = math.Sqrt(u*u+v*v) + C
		return
	}
	for i := 0; i < Nedge; i++ {
		iL := i + shiftL
		iR := Nedge - 1 - i + shiftR // Shared edges run in reverse order relative to each other
		qL := c.GetQQ(kL, KmaxL, iL, qfDL)
		rhoL, uL, vL = qL[0], qL[1]/qL[0], qL[2]/qL[0]
		pL, CL = c.FS.GetFlowFunction(qL, StaticPressure), c.FS.GetFlowFunction(qL, SoundSpeed)
		qR := c.GetQQ(kR, KmaxR, iR, qfDR)
		rhoR, uR, vR = qR[0], qR[1]/qR[0], qR[2]/qR[0]
		pR, CR = c.FS.GetFlowFunction(qR, StaticPressure), c.FS.GetFlowFunction(qR, SoundSpeed)
		FxL, FyL := c.CalculateFlux(kL, KmaxL, iL, qfDL)
		FxR, FyR := c.CalculateFlux(kR, KmaxR, iR, qfDR) // Reverse the right edge to match
		maxV := math.Max(maxVF(uL, vL, pL, rhoL, CL), maxVF(uR, vR, pR, rhoR, CR))
		indL, indR := kL+iL*KmaxL, kR+iR*KmaxR
		for n := 0; n < 4; n++ {
			normalFlux[i][n] = 0.5 * (normal[0]*(FxL[n]+FxR[n]) + normal[1]*(FyL[n]+FyR[n]))
			normalFlux[i][n] += 0.5 * maxV * (qfDL[n][indL] - qfDR[n][indR])
			normalFluxReversed[Nedge-1-i][n] = -normalFlux[i][n]
		}
	}
}

func (c *Euler) RoeFlux(kL, kR, KmaxL, KmaxR, shiftL, shiftR int,
	Q_FaceL, Q_FaceR [4]utils.Matrix, normal [2]float64, normalFlux, normalFluxReversed [][4]float64) {
	var (
		Nedge            = c.dfr.FluxElement.Nedge
		rhoL, uL, vL, pL float64
		rhoR, uR, vR, pR float64
		hL, hR           float64
		Gamma            = c.FS.Gamma
		GM1              = Gamma - 1
		qfDL, qfDR       = Get4DP(Q_FaceL), Get4DP(Q_FaceR)
	)
	rotateMomentum := func(k, Kmax, i int, qfD [4][]float64) {
		ind := k + i*Kmax
		um, vm := qfD[1][ind], qfD[2][ind]
		qfD[1][ind] = um*normal[0] + vm*normal[1]
		qfD[2][ind] = -um*normal[1] + vm*normal[0]
	}
	for i := 0; i < Nedge; i++ {
		iL := i + shiftL
		iR := Nedge - 1 - i + shiftR // Shared edges run in reverse order relative to each other
		// Rotate the momentum into face normal coordinates before calculating fluxes
		rotateMomentum(kL, KmaxL, iL, qfDL)
		rotateMomentum(kR, KmaxR, iR, qfDR)
		qL := c.GetQQ(kL, KmaxL, iL, qfDL)
		rhoL, uL, vL = qL[0], qL[1]/qL[0], qL[2]/qL[0]
		pL = c.FS.GetFlowFunction(qL, StaticPressure)
		qR := c.GetQQ(kR, KmaxR, iR, qfDR)
		rhoR, uR, vR = qR[0], qR[1]/qR[0], qR[2]/qR[0]
		pR = c.FS.GetFlowFunction(qR, StaticPressure)
		FxL, _ := c.CalculateFlux(kL, KmaxL, iL, qfDL)
		FxR, _ := c.CalculateFlux(kR, KmaxR, iR, qfDR) // Reverse the right edge to match
		/*
		   HM = (EnerM+pM).dd(rhoM);  HP = (EnerP+pP).dd(rhoP);
		*/
		// Enthalpy
		hL, hR = c.FS.GetFlowFunction(qL, Enthalpy), c.FS.GetFlowFunction(qR, Enthalpy)
		/*
			rhoMs = sqrt(rhoM); rhoPs = sqrt(rhoP);
			rhoMsPs = rhoMs + rhoPs;

			rho = rhoMs.dm(rhoPs);
			u   = (rhoMs.dm(uM) + rhoPs.dm(uP)).dd(rhoMsPs);
			v   = (rhoMs.dm(vM) + rhoPs.dm(vP)).dd(rhoMsPs);
			H   = (rhoMs.dm(HM) + rhoPs.dm(HP)).dd(rhoMsPs);
			c2 = gm1 * (H - 0.5*(sqr(u)+sqr(v)));
			c = sqrt(c2);
		*/
		// Compute Roe average variables
		rhoLs, rhoRs := math.Sqrt(rhoL), math.Sqrt(rhoR)
		rhoLsRs := rhoLs + rhoRs

		rho := rhoLs * rhoRs
		u := (rhoLs*uL + rhoRs*uR) / rhoLsRs
		v := (rhoLs*vL + rhoRs*vR) / rhoLsRs
		h := (rhoLs*hL + rhoRs*hR) / rhoLsRs
		c2 := GM1 * (h - 0.5*(u*u+v*v))
		c := math.Sqrt(c2)
		/*
		   dW1 = -0.5*(rho.dm(uP-uM)).dd(c) + 0.5*(pP-pM).dd(c2);
		   dW2 = (rhoP-rhoM) - (pP-pM).dd(c2);
		   dW3 = rho.dm(vP-vM);
		   dW4 = 0.5*(rho.dm(uP-uM)).dd(c) + 0.5*(pP-pM).dd(c2);

		   dW1 = abs(u-c).dm(dW1);
		   dW2 = abs(u  ).dm(dW2);
		   dW3 = abs(u  ).dm(dW3);
		   dW4 = abs(u+c).dm(dW4);
		*/
		// Riemann fluxes
		dW1 := -0.5*(rho*(uR-uL))/c + 0.5*(pR-pL)/c2
		dW2 := (rhoR - rhoL) - (pR-pL)/c2
		dW3 := rho * (vR - vL)
		dW4 := 0.5*(rho*(uR-uL))/c + 0.5*(pR-pL)/c2
		dW1 = math.Abs(u-c) * dW1
		dW2 = math.Abs(u) * dW2
		dW3 = math.Abs(u) * dW3
		dW4 = math.Abs(u+c) * dW4
		/*
		   DMat fx = (fxQP+fxQM)/2.0;
		   fx(All,1) -= (dW1               + dW2                                   + dW4              )/2.0;
		   fx(All,2) -= (dW1.dm(u-c)       + dW2.dm(u)                             + dW4.dm(u+c)      )/2.0;
		   fx(All,3) -= (dW1.dm(v)         + dW2.dm(v)                 + dW3       + dW4.dm(v)        )/2.0;
		   fx(All,4) -= (dW1.dm(H-u.dm(c)) + dW2.dm(sqr(u)+sqr(v))/2.0 + dW3.dm(v) + dW4.dm(H+u.dm(c)))/2.0;
		*/
		// Form Roe Fluxes
		for n := 0; n < 4; n++ {
			normalFlux[i][n] = 0.5 * (FxL[n] + FxR[n]) // Ave of normal component of flux
		}
		normalFlux[i][0] -= 0.5 * (dW1 + dW2 + dW4)
		normalFlux[i][1] -= 0.5 * (dW1*(u-c) + dW2*u + dW4*(u+c))
		normalFlux[i][2] -= 0.5 * (dW1*v + dW2*v + dW3 + dW4*v)
		normalFlux[i][3] -= 0.5 * (dW1*(h-u*c) + 0.5*dW2*(u*u+v*v) + dW3*v + dW4*(h+u*c))
		/*
		   flux = fx;    fx2.borrow(Ngf, fx.pCol(2)); fx3.borrow(Ngf, fx.pCol(3));
		   flux(All,2) = lnx.dm(fx2) - lny.dm(fx3);
		   flux(All,3) = lny.dm(fx2) + lnx.dm(fx3);
		*/
		// rotate back to Cartesian
		normalFlux[i][1], normalFlux[i][2] = normal[0]*normalFlux[i][1]-normal[1]*normalFlux[i][2], normal[1]*normalFlux[i][1]+normal[0]*normalFlux[i][2]
		for n := 0; n < 4; n++ {
			normalFluxReversed[Nedge-1-i][n] = -normalFlux[i][n]
		}
	}
}

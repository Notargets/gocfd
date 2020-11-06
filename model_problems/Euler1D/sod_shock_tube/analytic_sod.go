package sod_shock_tube

import (
	"math"
)

type State struct {
	rho, p, u, gamma float64
}

func (s State) C() float64 {
	return math.Sqrt(s.gamma * s.p / s.rho)
}

type SOD_Exact struct {
	gamma                            float64
	x_min, x_max, x0, x1, x2, x3, x4 float64
	rho_middle                       float64
	t                                float64
	l_s, r_s, post_s                 State
}

func NewSOD(t float64) (se *SOD_Exact) {
	gamma := 1.4
	se = &SOD_Exact{
		gamma: gamma,
		x_min: 0.,
		x_max: 1.,
		l_s:   State{1, 1, 0, gamma},
		r_s:   State{0.125, 0.1, 0, gamma},
	}
	se.Calc(t)
	return
}

func (se *SOD_Exact) Calc(t float64) {
	var (
		gamma    = se.gamma
		mu       = math.Sqrt((gamma - 1) / (gamma + 1))
		mu2      = mu * mu
		l_s, r_s = se.l_s, se.r_s
		P_post   = fzero(sod_func, math.Pi)
	)
	se.t = t
	se.rho_middle = l_s.rho * math.Pow(P_post/l_s.p, 1./gamma)
	se.x0 = 0.5 * (se.x_max + se.x_min)
	se.post_s = State{
		r_s.rho * ((P_post / r_s.p) + mu2) / (1 + mu2*(P_post/r_s.p)),
		P_post,
		r_s.u + (P_post-r_s.p)/math.Sqrt(0.5*r_s.rho*((gamma+1)*P_post+(gamma-1)*r_s.p)),
		gamma,
	}
	v_shock := se.post_s.u * (se.post_s.rho / r_s.rho) / ((se.post_s.rho / r_s.rho) - 1.)
	se.x1 = se.x0 - l_s.C()*t
	se.x3 = se.x0 + se.post_s.u*t
	se.x4 = se.x0 + v_shock*t
	//determining x2
	c_2 := l_s.C() - 0.5*(gamma-1.)*se.post_s.u
	se.x2 = se.x0 + t*(se.post_s.u-c_2)
}

func (se *SOD_Exact) Getx(x float64) (rho, p, u, e, rhou float64) {
	var (
		gamma              = se.gamma
		mu                 = math.Sqrt((gamma - 1) / (gamma + 1))
		mu2                = mu * mu
		l_s, r_s           = se.l_s, se.r_s
		x0, x1, x2, x3, x4 = se.x0, se.x1, se.x2, se.x3, se.x4
	)
	switch {
	case x < x1:
		rho = l_s.rho
		p = l_s.p
		u = l_s.u
	case x1 <= x && x <= x2:
		c := mu2*((x0-x)/se.t) + (1.-mu2)*l_s.C()
		rho = l_s.rho * math.Pow((c/l_s.C()), 2/(gamma-1))
		p = l_s.p * math.Pow(rho/l_s.rho, gamma)
		u = (1. - mu2) * ((-(x0 - x) / se.t) + l_s.C())
	case x2 <= x && x <= x3:
		rho = se.rho_middle
		p = se.post_s.p
		u = se.post_s.u
	case x3 <= x && x <= x4:
		rho = se.post_s.rho
		p = se.post_s.p
		u = se.post_s.u
	case x4 < x:
		rho = r_s.rho
		p = r_s.p
		u = r_s.u
	}
	e = p/(gamma-1.) + 0.5*u*u*rho
	rhou = rho * u
	return
}

func (se *SOD_Exact) Get() (X, Rho, P, RhoU, E []float64) {
	var (
		x_min, x_max, x1, x2, x3, x4 = se.x_min, se.x_max, se.x1, se.x2, se.x3, se.x4
	)
	tol := 0.0001
	midStep := (se.x2 - se.x1) / 10.
	X = []float64{
		x_min,
		x1 - tol, x1 + tol,
		x1 + midStep, x1 + 2*midStep, x1 + 3*midStep, x1 + 4*midStep, x1 + 5*midStep, x1 + 6*midStep, x1 + 7*midStep, x1 + 8*midStep, x1 + 9*midStep, x1 + 10*midStep - 2.*tol,
		x2 - tol, x2 + tol,
		x3 - tol, x3 + tol,
		x4 - tol, x4 + tol,
		x_max,
	}
	Rho = make([]float64, len(X))
	P = make([]float64, len(X))
	U := make([]float64, len(X))
	RhoU = make([]float64, len(X))
	E = make([]float64, len(X))
	for i, x := range X {
		Rho[i], P[i], U[i], E[i], RhoU[i] = se.Getx(x)
	}
	return
}

func fzero(f func(P float64) (y float64), start float64) float64 {
	var (
		tol = 0.0000001
		res float64
	)
	start_old := start / 2
	res = f(start_old)
	for math.Abs(res) > tol {
		resNew := f(start)
		deriv := (start - start_old) / (resNew - res)
		start_new := math.Abs(start - resNew*deriv)
		start_old = start
		start = start_new
		res = resNew
	}
	return start
}

func sod_func(P float64) (y float64) {
	var (
		rho_l, P_l = 1., 1.
		//rho_l, P_l, u_l = 1., 1., 0.
		//rho_r, P_r, u_r = 0.125, 0.1, 0.
		rho_r, P_r = 0.125, 0.1
		gamma      = 1.4
		mu         = math.Sqrt((gamma - 1) / (gamma + 1))
		mu2        = mu * mu
	)
	// Broken, per comments in Mathematica's community at https://www.mathworks.com/matlabcentral/fileexchange/46311-sod-shock-tube-problem-solver
	//y = (P-P_r)*math.Sqrt(utils.POW(1-mu2, 2)/(rho_r*(P+mu2*P_r))) - 2*(math.Sqrt(gamma)/(gamma-1))*(1-math.Pow(P, (gamma-1)/(2*gamma)))
	y = (P-P_r)*math.Sqrt((1-mu2)/(rho_r*(P+mu2*P_r))) - (math.Pow(P_l, (gamma-1)/(2*gamma))-math.Pow(P, (gamma-1)/(2*gamma)))*math.Sqrt(((1-mu2*mu2)*math.Pow(P_l, 1/gamma))/(mu2*mu2*rho_l))
	return
}

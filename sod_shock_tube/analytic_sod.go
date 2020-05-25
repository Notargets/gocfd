package sod_shock_tube

import (
	"fmt"
	"math"

	"github.com/notargets/gocfd/utils"
)

func SOD_calc(t float64) (X, Rho, P, U, E []float64) {
	var (
		x_min, x_max        = 0., 1.
		x0, rho_l, P_l, u_l = 0.5 * (x_max + x_min), 1., 1., 0.
		rho_r, P_r, u_r     = 0.125, 0.1, 0.
		gamma               = 1.4
		mu                  = math.Sqrt((gamma - 1) / (gamma + 1))
		c_l, c_r            = math.Sqrt(gamma * P_l / rho_l), math.Sqrt(gamma * P_r / rho_r)
		P_post              = fzero(sod_func, math.Pi)
		v_post              = 2 * (math.Sqrt(gamma) / (gamma - 1)) * (1 - math.Pow(P_post, (gamma-1)/(2*gamma)))
		rho_post            = rho_r * (((P_post / P_r) + mu*mu) / (1 + mu*mu*(P_post/P_r)))
		v_shock             = v_post * (rho_post / rho_r) / ((rho_post / rho_r) - 1.)
		rho_middle          = rho_l * math.Pow(P_post/P_l, 1./gamma)
		//Key Positions
		x1 = x0 - c_l*t
		x3 = x0 + v_post*t
		x4 = x0 + v_shock*t
		//determining x2
		c_2 = c_l - 0.5*(gamma-1.)*v_post
		x2  = x0 + t*(v_post-c_2)
		//boundaries (can be set)
	)
	fmt.Printf("Sod P_post = %v, sod_func(P_post) = %v\n", P_post, sod_func(P_post))
	_ = c_r
	tol := 0.00000001
	X = []float64{
		x_min,
		x1 - tol, x1 + tol,
		x2 - tol, x2 + tol,
		x3 - tol, x3 + tol,
		x4 - tol, x4 + tol,
		x_max,
	}
	Rho = make([]float64, len(X))
	P = make([]float64, len(X))
	U = make([]float64, len(X))
	E = make([]float64, len(X))
	for i, x := range X {
		switch {
		case x < x1:
			Rho[i] = rho_l
			P[i] = P_l
			U[i] = u_l
		case x1 <= x && x <= x2:
			c := mu*mu*((x0-X[i])/t) + (1.-mu*mu)*c_l
			Rho[i] = rho_l * math.Pow((c/c_l), 2/(gamma-1))
			P[i] = P_l * math.Pow(Rho[i]/rho_l, gamma)
			U[i] = (1. - mu*mu) * ((-(x0 - X[i]) / t) + c_l)
		case x2 <= x && x <= x3:
			Rho[i] = rho_middle
			P[i] = P_post
			U[i] = v_post
		case x3 <= x && x <= x4:
			Rho[i] = rho_post
			P[i] = P_post
			U[i] = v_post
		case x4 < x:
			Rho[i] = rho_r
			P[i] = P_r
			U[i] = u_r
		}
		E[i] = P[i] / ((gamma - 1.) * Rho[i])
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
		start_new := math.Abs(start - 0.01*f(start)/deriv)
		start_old = start
		start = start_new
		res = resNew
	}
	return start
}

func sod_func(P float64) (y float64) {
	var (
		//		rho_l, P_l, u_l = 1., 1., 0.
		//		rho_r, P_r, u_r = 0.125, 0.1, 0.
		rho_r, P_r = 0.125, 0.1
		gamma      = 1.4
		mu         = math.Sqrt((gamma - 1) / (gamma + 1))
		mu2        = mu * mu
	)
	y = (P-P_r)*math.Sqrt(utils.POW(1-mu2, 2)/(rho_r*(P+mu2*P_r))) - 2*(math.Sqrt(gamma)/(gamma-1))*(1-math.Pow(P, (gamma-1)/(2*gamma)))
	return
}

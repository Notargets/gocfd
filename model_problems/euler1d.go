package model_problems

import (
	"fmt"
	"math"
	"sync"
	"time"

	"github.com/notargets/avs/chart2d"
	utils2 "github.com/notargets/avs/utils"

	"github.com/notargets/gocfd/DG1D"
	"github.com/notargets/gocfd/utils"
)

type Euler1D struct {
	// Input parameters
	CFL, FinalTime  float64
	El              *DG1D.Elements1D
	RHSOnce         sync.Once
	Rho, RhoU, Ener utils.Matrix
}

func NewEuler1D(CFL, FinalTime float64, N, K int) (c *Euler1D) {
	VX, EToV := DG1D.SimpleMesh1D(-2, 2, K)
	c = &Euler1D{
		CFL:       CFL,
		FinalTime: FinalTime,
		El:        DG1D.NewElements1D(N, VX, EToV),
	}
	return
}

func (c *Euler1D) Run(showGraph bool, graphDelay ...time.Duration) {
	var (
		chart                  *chart2d.Chart2D
		colorMap               *utils2.ColorMap
		chartNameE, chartNameH string
		el                     = c.El
		resRho                 = utils.NewMatrix(el.Np, el.K)
		resRhoU                = utils.NewMatrix(el.Np, el.K)
		resEner                = utils.NewMatrix(el.Np, el.K)
		logFrequency           = 50
	)
	_, _, _ = resRho, resRhoU, resEner
	if showGraph {
		chart = chart2d.NewChart2D(1920, 1280, float32(el.X.Min()), float32(el.X.Max()), -1, 1)
		colorMap = utils2.NewColorMap(-1, 1, 1)
		chartNameE = "E Field"
		chartNameH = "H Field"
		go chart.Plot()
	}
	xmin := el.X.Row(1).Subtract(el.X.Row(0)).Apply(math.Abs).Min()
	dt := xmin * c.CFL
	Nsteps := int(math.Ceil(c.FinalTime / dt))
	dt = c.FinalTime / float64(Nsteps)
	fmt.Printf("FinalTime = %8.4f, Nsteps = %d, dt = %8.6f\n", c.FinalTime, Nsteps, dt)

	var Time float64
	for tstep := 0; tstep < Nsteps; tstep++ {
		if showGraph {
			if err := chart.AddSeries(chartNameE, el.X.Transpose().RawMatrix().Data, c.Rho.Transpose().RawMatrix().Data,
				chart2d.CrossGlyph, chart2d.Dashed, colorMap.GetRGB(0)); err != nil {
				panic("unable to add graph series")
			}
			if err := chart.AddSeries(chartNameH, el.X.Transpose().RawMatrix().Data, c.RhoU.Transpose().RawMatrix().Data,
				chart2d.CrossGlyph, chart2d.Dashed, colorMap.GetRGB(0.7)); err != nil {
				panic("unable to add graph series")
			}
			if err := chart.AddSeries(chartNameH, el.X.Transpose().RawMatrix().Data, c.Ener.Transpose().RawMatrix().Data,
				chart2d.CrossGlyph, chart2d.Dashed, colorMap.GetRGB(-0.7)); err != nil {
				panic("unable to add graph series")
			}
			if len(graphDelay) != 0 {
				time.Sleep(graphDelay[0])
			}
		}
		for INTRK := 0; INTRK < 5; INTRK++ {
			rhsRho, rhsRhoU, rhsEner := c.RHS()
			_, _, _ = rhsRho, rhsRhoU, rhsEner
		}
		Time += dt
		if tstep%logFrequency == 0 {
			fmt.Printf("Time = %8.4f, max_resid[%d] = %8.4f, emin = %8.6f, emax = %8.6f\n", Time, tstep, resEner.Max(), c.Ener.Min(), c.Ener.Max())
		}
	}
	return
}

func (c *Euler1D) RHS() (rhsRho, rhsRhoU, rhsEner utils.Matrix) {
	var (
		el       = c.El
		nrF, ncF = el.Nfp * el.NFaces, el.K
		gamma    = 1.4
		U        = c.RhoU.Copy().ElDiv(c.Rho)                     // Velocity
		RhoU2    = c.RhoU.Copy().POW(2).Scale(0.5).ElDiv(c.Rho)   // 1/2 * Rho * U^2
		Pres     = c.Ener.Copy().Subtract(RhoU2).Scale(gamma - 1) // Pressure
		CVel     = Pres.Copy().ElDiv(c.Rho).Scale(gamma).Apply(math.Sqrt)
		LM       = U.Copy().Apply(math.Abs).Add(CVel)
		RhoF     = c.RhoU.Copy()
		RhoUF    = RhoU2.Scale(2).Add(Pres)
		EnerF    = c.Ener.Copy().Add(Pres).ElDiv(c.RhoU.Copy().ElDiv(c.Rho))
		// Face jumps in primary and flux variables
		dRho   = c.Rho.Subset(el.VmapM, nrF, ncF).Subtract(c.Rho.Subset(el.VmapP, nrF, ncF))
		dRhoU  = c.RhoU.Subset(el.VmapM, nrF, ncF).Subtract(c.RhoU.Subset(el.VmapP, nrF, ncF))
		dEner  = c.Ener.Subset(el.VmapM, nrF, ncF).Subtract(c.Ener.Subset(el.VmapP, nrF, ncF))
		dRhoF  = RhoF.Subset(el.VmapM, nrF, ncF).Subtract(RhoF.Subset(el.VmapP, nrF, ncF))
		dRhoUF = RhoUF.Subset(el.VmapM, nrF, ncF).Subtract(RhoUF.Subset(el.VmapP, nrF, ncF))
		dEnerF = EnerF.Subset(el.VmapM, nrF, ncF).Subtract(EnerF.Subset(el.VmapP, nrF, ncF))
	)
	_, _, _, _, _, _, _ = LM, dRho, dRhoU, dEner, dRhoF, dRhoUF, dEnerF
	// Homogeneous boundary conditions at the inflow faces, Ez = 0
	// Reflection BC - Metal boundary - E is zero at shell face, H passes through (Neumann)
	// E on the boundary face is negative of E inside, so the diff in E at the boundary face is 2E of the interior
	//dE.AssignVector(el.MapB, c.E.SubsetVector(el.VmapB).Scale(2))
	// H on the boundary face is equal to H inside, so the diff in H at the boundary face is 0
	//dH.AssignVector(el.MapB, c.H.SubsetVector(el.VmapB).Set(0))

	// Upwind fluxes
	//fluxE = c.ZimPM.Copy().Add(c.ZimPP).POW(-1).ElMul(el.NX.Copy().ElMul(c.ZimPP).ElMul(dH).Subtract(dE))
	//fluxH = c.YimPM.Copy().Add(c.YimPP).POW(-1).ElMul(el.NX.Copy().ElMul(c.YimPP).ElMul(dE).Subtract(dH))

	//RHSE = el.Rx.Copy().Scale(-1).ElMul(el.Dr.Mul(c.H)).Add(el.LIFT.Mul(fluxE.ElMul(el.FScale))).ElDiv(c.Epsilon)
	//RHSH = el.Rx.Copy().Scale(-1).ElMul(el.Dr.Mul(c.E)).Add(el.LIFT.Mul(fluxH.ElMul(el.FScale))).ElDiv(c.Mu)

	return
}

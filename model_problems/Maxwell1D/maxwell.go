package Maxwell1D

import (
	"fmt"
	"math"
	"sync"
	"time"

	"github.com/notargets/gocfd/DG1D"
	"github.com/notargets/gocfd/utils"
)

type Maxwell struct {
	// Input parameters
	CFL, FinalTime                   float64
	El                               *DG1D.Elements1D
	RHSOnce, PlotOnce                sync.Once
	E, H                             utils.Matrix
	Epsilon, Mu                      utils.Matrix
	Zimp, ZimPM, ZimPP, YimPM, YimPP utils.Matrix
	ZimpDenom, YimpDenom             utils.Matrix
	model                            ModelType
}

type ModelType uint

const (
	GK ModelType = iota
	DFR
)

var (
	model_names = []string{
		"Galerkin Integration, Lax Flux",
		"DFR Integration, Lax Flux",
	}
)

func NewMaxwell(CFL, FinalTime float64, N, K int, model ModelType) (c *Maxwell) {
	VX, EToV := DG1D.SimpleMesh1D(-2, 2, K)
	c = &Maxwell{
		CFL:       CFL,
		FinalTime: FinalTime,
		El:        DG1D.NewElements1D(N, VX, EToV),
		model:     model,
	}
	fmt.Printf("CFL = %8.4f, Polynomial Degree N = %d (1 is linear), Num Elements K = %d\nModel Type: %s\n\n", CFL, N, K, model_names[c.model])
	epsData := utils.ConstArray(c.El.K, 1)
	ones := utils.NewVectorConstant(c.El.Np, 1)
	for i := c.El.K / 2; i < c.El.K; i++ {
		epsData[i] = 2
	}
	Eps1 := utils.NewVector(c.El.K, epsData)
	c.Epsilon = Eps1.Outer(ones)
	Mu1 := utils.NewVectorConstant(c.El.K, 1)
	c.Mu = Mu1.Outer(ones)
	c.E = c.El.X.Copy().Apply(func(val float64) float64 {
		if val < 0 {
			return math.Sin(math.Pi * val)
		} else {
			return 0
		}
	})
	c.H = utils.NewMatrix(c.El.Np, c.El.K)
	c.Zimp = c.Epsilon.Copy().POW(-1).ElMul(c.Mu).Apply(math.Sqrt)
	nrF, ncF := c.El.Nfp*c.El.NFaces, c.El.K
	c.ZimPM = c.Zimp.Subset(c.El.VmapM, nrF, ncF)
	c.ZimPP = c.Zimp.Subset(c.El.VmapP, nrF, ncF)
	c.ZimPM.SetReadOnly("ZimPM")
	c.ZimPP.SetReadOnly("ZimPP")
	c.YimPM, c.YimPP = c.ZimPM.Copy().POW(-1), c.ZimPP.Copy().POW(-1)
	c.YimPM.SetReadOnly("YimPM")
	c.YimPP.SetReadOnly("YimPP")
	return
}

func (c *Maxwell) Run(showGraph bool, graphDelay ...time.Duration) {
	var (
		el           = c.El
		resE         = utils.NewMatrix(el.Np, el.K)
		resH         = utils.NewMatrix(el.Np, el.K)
		logFrequency = 50
		rhs          func() (rhsE, rhsH utils.Matrix)
	)
	switch c.model {
	case GK:
		rhs = c.RHS_GK
	case DFR:
		fallthrough
	default:
		rhs = c.RHS_DFR
	}
	xmin := el.X.Row(1).Subtract(el.X.Row(0)).Apply(math.Abs).Min()
	dt := xmin * c.CFL
	Nsteps := int(math.Ceil(c.FinalTime / dt))
	dt = c.FinalTime / float64(Nsteps)
	fmt.Printf("FinalTime = %8.4f, Nsteps = %d, dt = %8.6f\n", c.FinalTime, Nsteps, dt)

	var Time float64
	for tstep := 0; tstep < Nsteps; tstep++ {
		for INTRK := 0; INTRK < 5; INTRK++ {
			rhsE, rhsH := rhs()
			resE.Scale(utils.RK4a[INTRK]).Add(rhsE.Scale(dt))
			resH.Scale(utils.RK4a[INTRK]).Add(rhsH.Scale(dt))
			c.E.Add(resE.Copy().Scale(utils.RK4b[INTRK]))
			c.H.Add(resH.Copy().Scale(utils.RK4b[INTRK]))
		}
		Time += dt
		if tstep%logFrequency == 0 {
			fmt.Printf("Time = %8.4f, max_resid[%d] = %8.4f, emin = %8.6f, emax = %8.6f\n", Time, tstep, resE.Max(), c.E.Min(), c.E.Max())
		}
	}
	return
}

func (c *Maxwell) RHS_DFR() (RHSE, RHSH utils.Matrix) {
	var (
		el       = c.El
		nrF, ncF = el.Nfp * el.NFaces, el.K
		// FluxE, FluxH = c.E.Copy(), c.H.Copy()
		FluxE, FluxH = c.E, c.H
		// Field flux differerence across faces
		dE                   = FluxE.Subset(el.VmapM, nrF, ncF).Subtract(FluxE.Subset(el.VmapP, nrF, ncF))
		dH                   = FluxH.Subset(el.VmapM, nrF, ncF).Subtract(FluxH.Subset(el.VmapP, nrF, ncF))
		FaceFluxE, FaceFluxH utils.Matrix
		aDiss2, aDiss4       = .03, 0.02
	)
	c.RHSOnce.Do(func() {
		c.ZimpDenom = c.ZimPM.Copy().Add(c.ZimPP).POW(-1)
		c.YimpDenom = c.YimPM.Copy().Add(c.YimPP).POW(-1)
	})
	// FluxE, FluxH = el.SlopeLimitN(FluxE, 20), el.SlopeLimitN(FluxH, 20)
	// Homogeneous boundary conditions at the inflow faces, Ez = 0
	// Reflection BC - Metal boundary - E is zero at shell face, H passes through (Neumann)
	// E on the boundary face is negative of E inside, so the diff in E at the boundary face is 2E of the interior
	dE.AssignVector(el.MapB, c.E.SubsetVector(el.VmapB).Scale(2))
	// H on the boundary face is equal to H inside, so the diff in H at the boundary face is 0
	dH.AssignVector(el.MapB, c.H.SubsetVector(el.VmapB).Set(0))

	fluxType := "anisotropic"
	switch fluxType {
	case "average":
		FaceFluxH = FluxH.Subset(el.VmapM, nrF, ncF).Add(FluxH.Subset(el.VmapP, nrF, ncF)).Scale(0.5)
		FaceFluxE = FluxE.Subset(el.VmapM, nrF, ncF).Add(FluxE.Subset(el.VmapP, nrF, ncF)).Scale(0.5)
	case "isotropic":
		FaceFluxH = FluxH.Subset(el.VmapM, nrF, ncF).Add(FluxH.Subset(el.VmapP, nrF, ncF)).Add(el.NX.Copy().ElMul(dE).ElMul(c.YimPP)).Scale(0.5)
		FaceFluxE = FluxE.Subset(el.VmapM, nrF, ncF).Add(FluxE.Subset(el.VmapP, nrF, ncF)).Add(el.NX.Copy().ElMul(dH).ElMul(c.ZimPP)).Scale(0.5)
	case "anisotropic":
		fallthrough
	default:
		FaceFluxH = FluxH.Subset(el.VmapM, nrF, ncF).ElMul(c.ZimPM).Add(FluxH.Subset(el.VmapP, nrF, ncF).ElMul(c.ZimPP)).Add(el.NX.Copy().ElMul(dE))
		FaceFluxH.ElMul(c.ZimpDenom)
		FaceFluxE = FluxE.Subset(el.VmapM, nrF, ncF).ElMul(c.YimPM).Add(FluxE.Subset(el.VmapP, nrF, ncF).ElMul(c.YimPP)).Add(el.NX.Copy().ElMul(dH))
		FaceFluxE.ElMul(c.YimpDenom)
	}

	FluxE.AssignVector(el.VmapM, FaceFluxE)
	FluxH.AssignVector(el.VmapM, FaceFluxH)

	GradE := el.Dr.Mul(FluxE)
	GradH := el.Dr.Mul(FluxH)
	GradE2 := el.Dr.Mul(GradE)
	GradH2 := el.Dr.Mul(GradH)
	var ADissE, ADissH utils.Matrix
	if true {
		// Ad-Hoc 2nd/4th Order Artificial Dissipation
		GradE4 := el.Dr.Mul(el.Dr.Mul(GradE2))
		GradH4 := el.Dr.Mul(el.Dr.Mul(GradH2))
		ADissE = GradE4.Scale(aDiss4).Add(GradE2.Scale(aDiss2))
		ADissH = GradH4.Scale(aDiss4).Add(GradH2.Scale(aDiss2))
	} else {
		ADissE = GradE2.Scale(aDiss2)
		ADissH = GradH2.Scale(aDiss2)
	}

	RHSE = GradH.ElDiv(c.Epsilon).Scale(-1).Add(ADissE).ElMul(el.Rx)
	RHSH = GradE.ElDiv(c.Mu).Scale(-1).Add(ADissH).ElMul(el.Rx)

	return
}

func (c *Maxwell) RHS_GK() (RHSE, RHSH utils.Matrix) {
	var (
		nrF, ncF = c.El.Nfp * c.El.NFaces, c.El.K
		// Field flux differerence across faces
		dE           = c.E.Subset(c.El.VmapM, nrF, ncF).Subtract(c.E.Subset(c.El.VmapP, nrF, ncF))
		dH           = c.H.Subset(c.El.VmapM, nrF, ncF).Subtract(c.H.Subset(c.El.VmapP, nrF, ncF))
		el           = c.El
		fluxE, fluxH utils.Matrix
	)
	// Homogeneous boundary conditions at the inflow faces, Ez = 0
	// Reflection BC - Metal boundary - E is zero at shell face, H passes through (Neumann)
	// E on the boundary face is negative of E inside, so the diff in E at the boundary face is 2E of the interior
	dE.AssignVector(el.MapB, c.E.SubsetVector(el.VmapB).Scale(2))
	// H on the boundary face is equal to H inside, so the diff in H at the boundary face is 0
	dH.AssignVector(el.MapB, c.H.SubsetVector(el.VmapB).Set(0))

	// Upwind fluxes
	fluxE = c.ZimPM.Copy().Add(c.ZimPP).POW(-1).ElMul(el.NX.Copy().ElMul(c.ZimPP).ElMul(dH).Subtract(dE))
	fluxH = c.YimPM.Copy().Add(c.YimPP).POW(-1).ElMul(el.NX.Copy().ElMul(c.YimPP).ElMul(dE).Subtract(dH))

	RHSE = el.Rx.Copy().Scale(-1).ElMul(el.Dr.Mul(c.H)).Add(el.LIFT.Mul(fluxE.ElMul(el.FScale))).ElDiv(c.Epsilon)
	RHSH = el.Rx.Copy().Scale(-1).ElMul(el.Dr.Mul(c.E)).Add(el.LIFT.Mul(fluxH.ElMul(el.FScale))).ElDiv(c.Mu)

	return
}

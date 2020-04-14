package model_problems

import (
	"math"
	"sync"
	"time"

	"github.com/notargets/gocfd/DG1D"
	"github.com/notargets/gocfd/utils"
)

type Maxwell1D struct {
	// Input parameters
	CFL, FinalTime                   float64
	El                               *DG1D.Elements1D
	RHSOnce                          sync.Once
	E, H                             utils.Matrix
	Epsilon, Mu                      utils.Matrix
	Zimp, ZimPM, ZimPP, YimPM, YimPP utils.Matrix
}

func NewMaxwell1D(CFL, FinalTime float64, N, K int) (c *Maxwell1D) {
	VX, EToV := DG1D.SimpleMesh1D(-1, 1, K)
	c = &Maxwell1D{
		CFL:       CFL,
		FinalTime: FinalTime,
		El:        DG1D.NewElements1D(N, VX, EToV),
	}
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
	c.Zimp = c.Epsilon.Copy().POW(-1).ElementMultiply(c.Mu).Apply(math.Sqrt)
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

func (c *Maxwell1D) Run(showGraph bool, graphDelay ...time.Duration) {
	var (
	//el        = c.El
	//chart     *chart2d.Chart2D
	//colorMap  *utils2.ColorMap
	//chartName string
	)
	c.RHS()
	return
}

func (c *Maxwell1D) RHS() (RHSE, RHSH utils.Matrix) {
	var (
		nrF, ncF = c.El.Nfp * c.El.NFaces, c.El.K
		// Field flux differerence across faces
		dE = c.E.Subset(c.El.VmapM, nrF, ncF).Subtract(c.E.Subset(c.El.VmapP, nrF, ncF))
		dH = c.H.Subset(c.El.VmapM, nrF, ncF).Subtract(c.H.Subset(c.El.VmapP, nrF, ncF))
		el = c.El
	)
	// Homogeneous boundary conditions, Ez = 0 (but note that this means dE is 2x the edge node value?)
	dE.AssignVector(c.El.MapB, c.E.SubsetVector(el.VmapB).Scale(2))
	dH.AssignVector(el.MapB, c.H.SubsetVector(el.VmapB).Set(0))

	// Upwind fluxes
	fluxE := c.ZimPM.Copy().Add(c.ZimPP).POW(-1).ElementMultiply(el.NX.Copy().ElementMultiply(c.ZimPP.Copy().ElementMultiply(dH).Subtract(dE)))
	fluxH := c.YimPM.Copy().Add(c.YimPP).POW(-1).ElementMultiply(el.NX.Copy().ElementMultiply(c.YimPP.Copy().ElementMultiply(dE).Subtract(dH)))

	RHSE = el.Rx.Copy().Scale(-1).ElementMultiply(el.Dr.Mul(c.H)).Add(el.LIFT.Mul(fluxE.ElementMultiply(el.FScale)).ElementDivide(c.Epsilon))
	RHSH = el.Rx.Copy().Scale(-1).ElementMultiply(el.Dr.Mul(c.E)).Add(el.LIFT.Mul(fluxH.ElementMultiply(el.FScale)).ElementDivide(c.Mu))

	return
}

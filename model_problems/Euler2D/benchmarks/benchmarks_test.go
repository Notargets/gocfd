package benchmarks

import (
	"fmt"
	"testing"

	"github.com/notargets/gocfd/model_problems/Euler2D"
)

func BenchmarkEulerSolve(b *testing.B) {
	var (
		plotMesh  = false
		pm        = &Euler2D.PlotMeta{Plot: false, StepsBeforePlot: 100}
		Nmax      = 2
		FinalTime = 0.1
		c         = make([]*Euler2D.Euler, Nmax+1)
	)
	for n := 1; n <= Nmax; n++ {
		c[n] = Euler2D.NewEuler(FinalTime, n, "../../../DG2D/vortex-new.su2", 1.00, Euler2D.FLUX_LaxFriedrichs, Euler2D.IVORTEX, 0, 0, 1.4, 0, false, 5000, plotMesh, false, false)
	}
	b.ResetTimer()
	// The benchmark loop
	for i := 0; i < b.N; i++ {
		// This is separate to enable easy performance and memory profiling
		for n := 1; n <= Nmax; n++ {
			c[n].Solve(pm)
		}
	}
}

func BenchmarkEulerGetFlowFunction(b *testing.B) {
	var (
		q   = [4]float64{1, 1, 1, 1}
		c   = Euler2D.NewEuler(1, 1, "", 1, Euler2D.FLUX_LaxFriedrichs, Euler2D.FREESTREAM, 1, 1, 1.4, 0, true, 1, false, false, false)
		GM1 = c.FS.Gamma - 1
	)
	var p float64
	b.Run("direct compute", func(b *testing.B) {
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			qq := 0.5 * (q[1]*q[1] + q[2]*q[2]) / q[0]
			p = GM1 * (q[3] + qq)
		}
	})
	b.Run("Optimized function call", func(b *testing.B) {
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			p = c.FS.GetFlowFunction(q, Euler2D.StaticPressure)
		}
	})
	pressFunc := func(q [4]float64) (p float64) {
		qq := 0.5 * (q[1]*q[1] + q[2]*q[2]) / q[0]
		p = GM1 * (q[3] + qq)
		return
	}
	b.Run("inline function call", func(b *testing.B) {
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			p = pressFunc(q)
		}
	})
	pressFunc2 := func(q0, q1, q2, q3 float64) (p float64) {
		qq := 0.5 * (q1*q1 + q2*q2) / q0
		p = GM1 * (q3 + qq)
		return
	}
	b.Run("inline function call, discrete args", func(b *testing.B) {
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			p = pressFunc2(q[0], q[1], q[2], q[3])
		}
	})

	fmt.Printf("p = %8.5f\n", p)
}

package DG2D

import (
	"fmt"
	"math"
	"strconv"
	"testing"
	"time"

	"github.com/notargets/avs/chart2d"
	utils2 "github.com/notargets/avs/utils"

	"github.com/notargets/gocfd/utils"
)

func _TestEdgeProjection(t *testing.T) {
	NMin := 1
	NMax := 1
	for N := NMin; N <= NMax; N++ {
		angle := 140.
		dfr := CreateEquiTriMesh(N, angle)
		NpInt := dfr.FluxElement.NpInt
		NpEdge := dfr.FluxElement.NpEdge
		// fmt.Printf("NInt[%d] = %d, NInt/3 = %d, Remainder=%d, NEdge=%d\n",
		// 	N, NpInt, NpInt/3, NpInt%3, NEdge)

		// The Edge Basis Vandermonde matrices are all diagonal,
		// as the edge basis is Lagrange. The diagonal element is 2
		tangentNeg := func(v [2]float64) (tt [2]float64) {
			// Left hand rule tangent to follow counter-clockwise edge traversal
			tt = [2]float64{-v[1], v[0]}
			tNorm := math.Sqrt(tt[0]*tt[0] + tt[1]*tt[1])
			tt[0] /= tNorm
			tt[1] /= tNorm
			return
		}
		rsFromTangentProj := func(r, s float64,
			normal, origin [2]float64) (rP, sP float64) {
			tt := tangentNeg(normal)
			dot := tt[0]*(r-origin[0]) + tt[1]*(s-origin[1])
			rP = tt[0]*dot + origin[0]
			sP = tt[1]*dot + origin[1]
			return
		}
		var AProj [3]utils.Matrix
		var A utils.Matrix
		for i := 0; i < 3; i++ {
			A = utils.NewMatrix(NpInt, NpEdge)
		}

		var origin [2]float64
		for nEdge := 0; nEdge < 3; nEdge++ {
			switch nEdge {
			case 0:
				origin = [2]float64{-1, -1}
			case 1:
				origin = [2]float64{1, -1}
			case 2:
				origin = [2]float64{-1, 1}
			}
			// fmt.Printf("nEdge: %d\n", nEdge)
			offset := 2*NpInt + nEdge*NpEdge
			for j := offset; j < offset+NpEdge; j++ {
				b_j := dfr.FluxElement.DOFVectors[j] // Edge normal vector
				// fmt.Printf("b_j: [%5.2f,%5.2f]\n", b_j.Eval()[0], b_j.Eval()[1])
				jj := j - offset
				for i := 0; i < NpInt; i++ {
					r, s := dfr.FluxElement.R.AtVec(i), dfr.FluxElement.S.AtVec(i)
					rP, sP := rsFromTangentProj(r, s, b_j.Eval(), origin)
					// fmt.Printf("Proj(%5.2f,%5.2f)=[%5.2f,%5.2f]\n", r, s, rP, sP)
					A.Set(i, jj, dfr.FluxElement.Phi[j].Dot(rP, sP, b_j.Eval()))
				}
			}
			A.Print("A")
			// B := A.Transpose().Mul(A)
			// B.Print("B")
			// B.InverseWithCheck().Print("BInv")
			AProj[nEdge] = A.Transpose().Mul(A).InverseWithCheck().Mul(A.Transpose())
			// AProj[nEdge].Print("Edge" + strconv.Itoa(nEdge))
		}
		_ = AProj
		U1, U2 := ShockConditions(1.1, 0)
		QSol, QFlux := SetShockConditions(dfr, U1, U2)
		QSol.Transpose().Print("QSol")
		_ = QFlux
		fmt.Printf("Number of edge points is: %d\n", NpEdge)
		Dens := QSol.Col(0)
		// for i := 0; i < NpInt; i++ {
		// 	Dens.Set(2.7)
		// }
		Dens.Transpose().Print("Dens")
		for i := 0; i < 3; i++ {
			// EdgeProj := AProj[i].Mul(Dens.ToMatrix()).Scale(0.5)
			EdgeProj := AProj[i].Mul(Dens.ToMatrix()).Scale(2.)
			EdgeProj.Transpose().Print("EdgeProj" + strconv.Itoa(i))
		}
		// QFlux.Print("QFlux")
		// Dens := CopyDensityFromQAndEdges(dfr, QSol, QFlux)
		// var leftRight [2][]int
		// for i := 0; i < NpInt; i++ {
		// 	x, y := dfr.SolutionX.At(i, 0), dfr.SolutionY.At(i, 0)
		// 	var label string
		// }
		// Dens.Print("Dens")
	}
}

func TestEdgeInterpolation(t *testing.T) {
	NMin := 7
	NMax := 7
	// for N := NMin; N <= NMax; N++ {
	for N := NMin; N <= NMax; N++ {
		angle := 140.
		dfr := CreateEquiTriMesh(N, angle)
		NpInt := dfr.FluxElement.NpInt
		// NEdge := dfr.FluxElement.NpEdge
		// fmt.Printf("NInt[%d] = %d, NInt/3 = %d, Remainder=%d, NEdge=%d\n",
		// 	N, NpInt, NpInt/3, NpInt%3, NEdge)
		M := GradientMassMatrix(dfr)
		_ = M
		U1, U2 := ShockConditions(2, 0)
		QSol, QFlux := SetShockConditions(dfr, U1, U2)
		// QSol.Print("QSol")
		// QFlux.Print("QFlux")
		MInv := M.Transpose().Mul(M).InverseWithCheck().Mul(M.Transpose())
		// MInv.Print("MInv")
		Dens := CopyDensityFromQAndEdges(dfr, QSol, QFlux)
		Grad := MInv.Mul(Dens)
		// Grad.Print("Grad")
		GradA, GradB, GradC := Grad.DataP[0], Grad.DataP[1], Grad.DataP[2]
		n1, n2 := GradB, GradC
		norm := math.Sqrt(n1*n1 + n2*n2)
		n1 /= norm
		n2 /= norm
		fmt.Printf("Order[%d] - Shock normal vector: [%5.2f,%5.2f]\n",
			N, n1, n2)
		// Now compute an estimate of the shock threshold rho_thresh as the
		// average of min/max density in the element
		rhoMin, rhoMax := Dens.Min(), Dens.Max()
		rho_thresh := (rhoMin + rhoMax) / 2.
		RhoGrad := func(x, y float64) (RhoFit float64) {
			// Quickly determines which side of the shock surface we are on
			RhoFit = GradA + GradB*x + GradC*y
			return
		}
		// If RhoGrad > rho_thresh, we're on the post shock side
		_, _ = rho_thresh, RhoGrad
		var leftRight [2][]int
		for i := 0; i < NpInt; i++ {
			x, y := dfr.SolutionX.At(i, 0), dfr.SolutionY.At(i, 0)
			rhoFit := RhoGrad(x, y)
			var label string
			if rhoFit > rho_thresh {
				label = "post"
				leftRight[1] = append(leftRight[1], i)
			} else {
				label = "pre"
				leftRight[0] = append(leftRight[0], i)
			}
			fmt.Printf("Solution Point[%d][%5.2f,"+
				"%5.2f] %s shock rho=[%5.2f], rhoFit=[%5.2f]\n",
				i, x, y, label, QSol.At(i, 0), rhoFit)
		}
		if testing.Verbose() {
			PlotShockPoints(dfr, leftRight, 30*time.Second)
		}
		// edgeGroups := groupInteriorPoints(dfr.FluxX.DataP, dfr.FluxY.DataP,
		// 	NpInt, NEdge)
		// fmt.Println(edgeGroups)
		// if testing.Verbose() {
		// 	PlotEdgeGroups(dfr, edgeGroups)
		// }
	}
}

func PlotShockPoints(dfr *DFR2D, leftRight [2][]int, delay time.Duration) {
	var Line []float32
	addVertex := func(x, y, width float64) {
		wo2 := float32(width / 2.)
		xl, yl := float32(x), float32(y)
		Line = append(Line, xl-wo2, yl, xl+wo2, yl)
		Line = append(Line, xl, yl-wo2, xl, yl+wo2)
	}
	addLineSegment := func(x1, y1, x2, y2 float64) {
		x1l, y1l := float32(x1), float32(y1)
		x2l, y2l := float32(x2), float32(y2)
		Line = append(Line, x1l, y1l, x2l, y2l)
	}

	ch := chart2d.NewChart2D(-1, 1, -1, 1, 1024, 1024, utils2.WHITE,
		utils2.BLACK, 0.8)
	Line = append(Line, 0, -1, 0, 1)
	ch.AddLine(Line, utils2.WHITE)
	Line = []float32{}
	// Np := dfr.FluxElement.Np
	for k := 0; k < dfr.K; k++ {
		for ii, position := range leftRight {
			for _, pt := range position {
				addVertex(dfr.SolutionX.At(pt, k), dfr.SolutionY.At(pt, k), 0.025)
			}
			switch ii {
			case 0:
				ch.AddLine(Line, utils2.RED)
			case 1:
				ch.AddLine(Line, utils2.GREEN)
			}
			Line = []float32{}
		}
		// for i := 0; i < Np; i++ {
		// 	addVertex(dfr.FluxX.At(i, k), dfr.FluxY.At(i, k), 0.025)
		// }
	}
	for k := 0; k < dfr.K; k++ {
		verts := dfr.Tris.GetTriVerts(uint32(k))
		for ii := 0; ii < 3; ii++ {
			x1, y1 := dfr.VX.DataP[verts[ii]], dfr.VY.DataP[verts[ii]]
			ip := ii + 1
			if ii == 2 {
				ip = 0
			}
			x2, y2 := dfr.VX.DataP[verts[ip]], dfr.VY.DataP[verts[ip]]
			addLineSegment(x1, y1, x2, y2)
			ch.AddLine(Line, utils2.WHITE)
			Line = []float32{}
		}
	}
	// ch.NewWindow("Unit Triangle", 0.9, screen.AUTO)
	// Line = []float32{}
	// addLineSegment(-1, -1, 1, -1)
	// addLineSegment(1, -1, -1, 1)
	// addLineSegment(-1, 1, -1, -1)
	// ch.AddLine(Line, utils2.WHITE)
	time.Sleep(delay)
}

func CopyDensityFromQAndEdges(dfr *DFR2D, QSol, QFlux utils.Matrix) (Dens utils.Matrix) {
	var (
		// Np     = dfr.SolutionElement.Np
		// NpFlux = dfr.FluxElement.Np
		NpInt  = dfr.FluxElement.NpInt
		NpEdge = dfr.FluxElement.NpEdge
		K      = 1
	)
	Dens = utils.NewMatrix(NpInt+3*NpEdge, K)
	for i := 0; i < NpInt; i++ {
		Dens.Set(i, 0, QSol.At(i, 0))
	}
	for i := NpInt; i < NpInt+3*NpEdge; i++ {
		Dens.Set(i, 0, QFlux.At(i+NpInt, 0))
	}
	return
}

func SetShockConditions(dfr *DFR2D, U1, U2 [4]float64) (QSol, QFlux utils.Matrix) {
	var (
		Np     = dfr.SolutionElement.Np
		NpFlux = dfr.FluxElement.Np
		// NpInt  = dfr.FluxElement.NpInt
		// NpEdge = dfr.FluxElement.NpEdge
		K = 1
	)
	QSol = utils.NewMatrix(Np, 4*K)
	QFlux = utils.NewMatrix(NpFlux, 4*K)
	var x float64
	for i := 0; i < NpFlux; i++ {
		if i < Np {
			x, _ = dfr.SolutionX.At(i, 0), dfr.SolutionY.At(i, 0)
		} else {
			x, _ = dfr.FluxX.At(i, 0), dfr.FluxY.At(i, 0)
		}
		// fmt.Printf("x,y[%d] = [%5.2f,%5.2f]\n", i, x, y)
		for n := 0; n < 4; n++ {
			if x < 0 {
				if i < Np {
					QSol.DataP[n+i*4*K] = U1[n]
				}
				QFlux.DataP[n+i*4*K] = U1[n]
			} else {
				if i < Np {
					QSol.DataP[n+i*4*K] = U2[n]
				}
				QFlux.DataP[n+i*4*K] = U2[n]
			}
		}
	}
	return
}

func GradientMassMatrix(dfr *DFR2D) (M utils.Matrix) {
	var (
		// NpFlux = dfr.FluxElement.Np
		NpInt  = dfr.FluxElement.NpInt
		NpEdge = dfr.FluxElement.NpEdge
	)
	// Compute gradient using all points in flux element with density
	// First compose a mass matrix with all of the point locations
	M = utils.NewMatrix(NpInt+3*NpEdge, 3) // Points x (1:X:Y)
	for i := 0; i < NpInt; i++ {
		M.Set(i, 0, 1.)
		M.Set(i, 1, dfr.FluxX.DataP[i])
		M.Set(i, 2, dfr.FluxY.DataP[i])
	}
	for i := NpInt; i < NpInt+3*NpEdge; i++ {
		M.Set(i, 0, 1.)
		M.Set(i, 1, dfr.FluxX.DataP[i+NpInt])
		M.Set(i, 2, dfr.FluxY.DataP[i+NpInt])
	}
	// M.Print("M")
	return
}

// projectToEdge projects a 2D point (x, y) onto the edge defined by (x0, y0) to (x1, y1)
// using the edge’s tangent direction.
func projectToEdge(x, y, x0, y0, x1, y1 float64) float64 {
	// Compute the unit tangent vector of the edge.
	tanX := x1 - x0
	tanY := y1 - y0
	length := math.Hypot(tanX, tanY)
	tanX /= length
	tanY /= length

	// The 1D coordinate is the dot product of (x-x0, y-y0) with the tangent.
	return (x-x0)*tanX + (y-y0)*tanY
}

func PlotEdgeGroups(dfr *DFR2D, edgeGroups [3][]int) {
	var Line []float32
	addVertex := func(x, y, width float64) {
		wo2 := float32(width / 2.)
		xl, yl := float32(x), float32(y)
		Line = append(Line, xl-wo2, yl, xl+wo2, yl)
		Line = append(Line, xl, yl-wo2, xl, yl+wo2)
	}
	addLineSegment := func(x1, y1, x2, y2 float64) {
		x1l, y1l := float32(x1), float32(y1)
		x2l, y2l := float32(x2), float32(y2)
		Line = append(Line, x1l, y1l, x2l, y2l)
	}

	ch := chart2d.NewChart2D(-1, 1, -1, 1, 1024, 1024, utils2.WHITE,
		utils2.BLACK, 0.8)
	// Np := dfr.FluxElement.Np
	for k := 0; k < dfr.K; k++ {
		for ii, edge := range edgeGroups {
			for _, pt := range edge {
				addVertex(dfr.FluxX.At(pt, k), dfr.FluxY.At(pt, k), 0.025)
			}
			switch ii {
			case 0:
				ch.AddLine(Line, utils2.RED)
			case 1:
				ch.AddLine(Line, utils2.GREEN)
			case 2:
				ch.AddLine(Line, utils2.BLUE)
			}
			Line = []float32{}
		}
		// for i := 0; i < Np; i++ {
		// 	addVertex(dfr.FluxX.At(i, k), dfr.FluxY.At(i, k), 0.025)
		// }
	}
	for k := 0; k < dfr.K; k++ {
		verts := dfr.Tris.GetTriVerts(uint32(k))
		for ii := 0; ii < 3; ii++ {
			x1, y1 := dfr.VX.DataP[verts[ii]], dfr.VY.DataP[verts[ii]]
			ip := ii + 1
			if ii == 2 {
				ip = 0
			}
			x2, y2 := dfr.VX.DataP[verts[ip]], dfr.VY.DataP[verts[ip]]
			addLineSegment(x1, y1, x2, y2)
			switch ii {
			case 0:
				ch.AddLine(Line, utils2.RED)
			case 1:
				ch.AddLine(Line, utils2.GREEN)
			case 2:
				ch.AddLine(Line, utils2.BLUE)
			}
			Line = []float32{}
		}
	}
	// ch.NewWindow("Unit Triangle", 0.9, screen.AUTO)
	// Line = []float32{}
	// addLineSegment(-1, -1, 1, -1)
	// addLineSegment(1, -1, -1, 1)
	// addLineSegment(-1, 1, -1, -1)
	// ch.AddLine(Line, utils2.WHITE)
	time.Sleep(30 * time.Second)
}

// candidate holds per‐interior point distance info.
type candidate struct {
	index     int
	distances [3]float64 // min distance to each edge group.
	primary   int        // index of closest edge.
	secondary int        // index of next closest edge.
	delta     float64    // difference between secondary and primary distance.
}

// groupInteriorPoints first assigns each interior point (indices 0..NpInt-1)
// to the edge with the minimum distance, then rebalances if any group exceeds its target.
// The edges are at:
//
//	Edge 1: indices [2*NpInt, 2*NpInt+NpEdge)
//	Edge 2: indices [2*NpInt+NpEdge, 2*NpInt+2*NpEdge)
//	Edge 3: indices [2*NpInt+2*NpEdge, 2*NpInt+3*NpEdge)
func groupInteriorPoints(x, y []float64, NpInt, NpEdge int) [3][]int {
	cands := make([]candidate, NpInt)
	// squaredDistance returns the squared Euclidean distance.
	squaredDistance := func(x1, y1, x2, y2 float64) float64 {
		dx := x1 - x2
		dy := y1 - y2
		return dx*dx + dy*dy
	}
	// For each interior point, compute its distance to each edge.
	for i := 0; i < NpInt; i++ {
		best := math.MaxFloat64
		second := math.MaxFloat64
		bestEdge := -1
		secondEdge := -1

		// For each of the three edges.
		for e := 0; e < 3; e++ {
			start := 2*NpInt + e*NpEdge
			end := start + NpEdge
			distE := math.MaxFloat64
			for j := start; j < end; j++ {
				d := squaredDistance(x[i], y[i], x[j], y[j])
				if d < distE {
					distE = d
				}
			}
			cands[i].distances[e] = distE
			// Determine best and second best edges.
			if distE < best {
				second = best
				secondEdge = bestEdge
				best = distE
				bestEdge = e
			} else if distE < second {
				second = distE
				secondEdge = e
			}
		}
		cands[i].index = i
		cands[i].primary = bestEdge
		cands[i].secondary = secondEdge
		cands[i].delta = second - best
	}

	// Initial assignment: each interior point goes to its primary edge.
	assigned := make([]int, NpInt)
	for i := 0; i < NpInt; i++ {
		assigned[i] = cands[i].primary
	}

	// Count assignments per edge.
	counts := [3]int{}
	for i := 0; i < NpInt; i++ {
		counts[assigned[i]]++
	}

	// Determine target counts.
	base := NpInt / 3
	rem := NpInt % 3
	targets := [3]int{}
	for e := 0; e < 3; e++ {
		if e < rem {
			targets[e] = base + 1
		} else {
			targets[e] = base
		}
	}

	// Rebalancing: for any point in an over-target group, if its secondary candidate is under-target,
	// reassign it. We do a simple greedy pass.
	for {
		moved := false
		for i := 0; i < NpInt; i++ {
			curGroup := assigned[i]
			// Only consider points in groups that have too many points.
			if counts[curGroup] > targets[curGroup] {
				sec := cands[i].secondary
				// If the secondary candidate exists and is under target, reassign.
				if sec >= 0 && counts[sec] < targets[sec] {
					assigned[i] = sec
					counts[curGroup]--
					counts[sec]++
					moved = true
					// Break out to restart the loop after a move.
					break
				}
			}
		}
		if !moved {
			break
		}
	}

	// Build final groups.
	var groups [3][]int
	for i, grp := range assigned {
		groups[grp] = append(groups[grp], i)
	}
	return groups
}

func TestPlotEquiTri(t *testing.T) {
	// dfr := CreateEquiTriMesh(1, 0.15)
	// dfr := CreateEquiTriMesh(2, 0.15)
	dfr := CreateEquiTriMesh(2, 0.15)
	_ = dfr
	// if testing.Verbose() {
	// 	PlotDFRElements(dfr)
	// }
	Mach := 1.8
	Alpha := 10.
	fmt.Printf("Testing Mach %5.2f, Shock Angle:%5.2f\n", Mach, Alpha)
	U1, U2 := ShockConditions(Mach, Alpha)
	PrintU("U1", U1)
	PrintU("U2", U2)
	K := 1
	Np := dfr.SolutionElement.Np
	Q := utils.NewMatrix(Np, 4*K)
	for i := 0; i < Np; i++ {
		x, _ := dfr.SolutionX.At(i, 0), dfr.SolutionY.At(i, 0)
		// fmt.Printf("x,y[%d] = [%5.2f,%5.2f]\n", i, x, y)
		for n := 0; n < 4; n++ {
			if x < 0 {
				Q.DataP[n+i*4*K] = U1[n]
			} else {
				Q.DataP[n+i*4*K] = U2[n]
			}
		}
	}
	// Q.Print("Q1")
	// dfr.FluxEdgeInterp.Print("FluxEdgeInterp")
	QEdge := dfr.FluxEdgeInterp.Mul(Q)
	// QEdge.Print("QEdge")
	// dfr.FluxElement.ProjectFunctionOntoDOF()
	NpFluxEdge := dfr.FluxElement.NpEdge
	var CorrectU [4]float64
	var RMS, Max, MaxPercent [4]float64
	for i := 0; i < NpFluxEdge*3; i++ {
		x, _ := dfr.FluxX.At(i, 0), dfr.FluxY.At(i, 0)
		// fmt.Printf("x,y[%d] = [%5.2f,%5.2f]\n", i, x, y)
		if x < 0 {
			CorrectU = U1
		} else {
			CorrectU = U2
		}
		for n := 0; n < 4; n++ {
			err := QEdge.At(i, n) - CorrectU[n]
			RMS[n] += err * err
			if math.Abs(err) > Max[n] {
				Max[n] = err
				if math.Abs(CorrectU[n]) < 0.00001 {
					MaxPercent[n] = 100. * err
				} else {
					MaxPercent[n] = 100. * err / CorrectU[n]
				}
			}
		}
	}
	for n := 0; n < 4; n++ {
		RMS[n] = math.Sqrt(RMS[n] / float64(NpFluxEdge*3))
		fmt.Printf("Error[%d]:RMS%+5.2f, Max:%+5.2f,%+5.1f%%\n",
			n, RMS[n], Max[n], MaxPercent[n])
	}
}

func PrintU(label string, U [4]float64) {
	fmt.Printf("%s:", label)
	for i := 0; i < 4; i++ {
		fmt.Printf("%0.2f ", U[i])
	}
	fmt.Printf("\n")
}

func TestMachConditions(t *testing.T) {
	U1, U2 := ShockConditions(2, 0)
	fmt.Printf("U1: %5.2f,%5.2f,%5.2f,%5.2f\n", U1[0], U1[1], U1[2], U1[3])
	fmt.Printf("U2: %5.2f,%5.2f,%5.2f,%5.2f\n", U2[0], U2[1], U2[2], U2[3])
	U1, U2 = ShockConditions(2, 15)
	fmt.Printf("U1: %5.2f,%5.2f,%5.2f,%5.2f\n", U1[0], U1[1], U1[2], U1[3])
	fmt.Printf("U2: %5.2f,%5.2f,%5.2f,%5.2f\n", U2[0], U2[1], U2[2], U2[3])
}

// ShockConditions computes post-shock properties given pre-shock Mach number (M1) and shock angle (alpha)
// where alpha = 0 corresponds to a normal shock.
func ShockConditions(M1, alpha float64) (U1, U2 [4]float64) {
	// Alpha = 0 for a normal shock.
	// Alpha can be used to determine oblique shock results in inviscid flow
	// accurately
	// For normal shock, set beta = alpha + 90 degrees.
	// Here, beta is the shock angle relative to the incoming flow.
	var (
		gamma   = 1.4
		beta    = alpha + 90.
		betaRad = beta * math.Pi / 180.0 // Convert beta to radians
	)

	// Compute normal component of Mach number.
	M1n := M1 * math.Sin(betaRad)

	// Compute density ratio from normal shock relations.
	rhoRatio := ((gamma + 1) * M1n * M1n) / ((gamma-1)*M1n*M1n + 2)

	// Compute pressure ratio from normal shock relations.
	pRatio := 1 + (2*gamma/(gamma+1))*(M1n*M1n-1)

	// Pre-shock conditions (assumed non-dimensional for testing).
	rho1 := 1.0 // Reference density.
	p1 := 1.0   // Reference pressure.
	u1 := M1    // Freestream velocity (assumed normalized so that a1 = 1).
	v1 := 0.0   // No initial vertical velocity.

	En1 := p1/(gamma-1) + 0.5*rho1*(u1*u1+v1*v1)
	U1 = [4]float64{rho1, u1, v1, En1}
	if M1 < 1. {
		U2 = U1
		return
	}
	// Define an inline Newton solver to compute u2n.
	NewtonSolver := func(rho1, u1n, p1, rho2, p2, tol float64, maxIter int) float64 {
		u2n := u1n / 2.0 // Initial guess.
		for i := 0; i < maxIter; i++ {
			// Correct residual: F = (rho2*u2n^2 + p2) - (rho1*u1n^2 + p1)
			F := (rho2*u2n*u2n + p2) - (rho1*u1n*u1n + p1)
			// Derivative: F' = 2 * rho2 * u2n.
			Fprime := 2.0 * rho2 * u2n

			// Newton update.
			u2nNew := u2n - F/Fprime

			// Check convergence.
			if math.Abs(u2nNew-u2n) < tol {
				return u2nNew
			}
			u2n = u2nNew
		}
		fmt.Println("Warning: Newton solver did not converge!")
		return u2n
	}

	// Compute post-shock properties.
	rho2 := rho1 * rhoRatio
	p2 := p1 * pRatio
	// Tangential component remains unchanged.
	u1t := u1 * math.Cos(betaRad)

	// Solve for post-shock normal velocity using the Newton solver.
	u2n := NewtonSolver(rho1, M1n, p1, rho2, p2, 1e-6, 50)

	// Compute full post-shock velocity components.
	u2 := u2n*math.Sin(betaRad) + u1t*math.Cos(betaRad)
	v2 := u2n*math.Cos(betaRad) - u1t*math.Sin(betaRad)

	// Compute total energy (using post-shock pressure and velocity).
	En2 := p2/(gamma-1) + 0.5*rho2*(u2*u2+v2*v2)

	// Pack pre- and post-shock states into U1 and U2.
	U2 = [4]float64{rho2, u2, v2, En2}
	return
}

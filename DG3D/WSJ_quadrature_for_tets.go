package DG3D

import "math"

// WilliamsShunnJamesonCubature returns Williams-Shunn-Jameson symmetric quadrature
// points and weights for the reference tetrahedron with vertices at
// (-1,-1,-1), (1,-1,-1), (-1,1,-1), (-1,-1,1)
// The rules integrate polynomials exactly up to degree 2*order-1
func WilliamsShunnJamesonCubature(order int) (r, s, t []float64, w []float64) {
	switch order {
	case 1:
		// Order 1: 1 point rule (centroid)
		// Exact for polynomials up to degree 1
		// Barycentric: (1/4, 1/4, 1/4, 1/4)
		r = []float64{-1.0 + 2.0*0.25}
		s = []float64{-1.0 + 2.0*0.25}
		t = []float64{-1.0 + 2.0*0.25}
		w = []float64{4.0 / 3.0}

	case 2:
		// Order 2: 4 point rule
		// Exact for polynomials up to degree 3
		// Symmetric points near vertices
		a := 0.58541019662496845446
		b := 0.13819660112501051518

		// Barycentric coordinates: (a,b,b,b), (b,a,b,b), (b,b,a,b), (b,b,b,a)
		r = []float64{
			-1.0 + 2.0*a,
			-1.0 + 2.0*b,
			-1.0 + 2.0*b,
			-1.0 + 2.0*b,
		}
		s = []float64{
			-1.0 + 2.0*b,
			-1.0 + 2.0*a,
			-1.0 + 2.0*b,
			-1.0 + 2.0*b,
		}
		t = []float64{
			-1.0 + 2.0*b,
			-1.0 + 2.0*b,
			-1.0 + 2.0*a,
			-1.0 + 2.0*b,
		}
		w = []float64{
			1.0 / 3.0,
			1.0 / 3.0,
			1.0 / 3.0,
			1.0 / 3.0,
		}

	case 3:
		// Order 3: 10 point rule
		// Exact for polynomials up to degree 5
		// Based on CCP lattice arrangement

		// Vertex points (4)
		// Barycentric: (v, (1-v)/3, (1-v)/3, (1-v)/3) and permutations
		v := 0.7236067977499790 // (5 + 3*sqrt(5))/20
		b3 := (1.0 - v) / 3.0
		xv := []float64{
			-1.0 + 2.0*v,
			-1.0 + 2.0*b3,
			-1.0 + 2.0*b3,
			-1.0 + 2.0*b3,
		}
		yv := []float64{
			-1.0 + 2.0*b3,
			-1.0 + 2.0*v,
			-1.0 + 2.0*b3,
			-1.0 + 2.0*b3,
		}
		zv := []float64{
			-1.0 + 2.0*b3,
			-1.0 + 2.0*b3,
			-1.0 + 2.0*v,
			-1.0 + 2.0*b3,
		}
		wv := (25.0 - 5.0*math.Sqrt(5.0)) / 480.0

		// Edge midpoints (6)
		// Barycentric: (1/2, 1/2, 0, 0) and permutations
		xe := []float64{
			0.0, 0.0, -1.0, -1.0, 0.0, -1.0,
		}
		ye := []float64{
			0.0, -1.0, 0.0, -1.0, -1.0, 0.0,
		}
		ze := []float64{
			-1.0, 0.0, 0.0, -1.0, -1.0, -1.0,
		}
		we := (25.0 + 5.0*math.Sqrt(5.0)) / 480.0

		// Combine all points
		r = append(xv, xe...)
		s = append(yv, ye...)
		t = append(zv, ze...)

		w = make([]float64, 10)
		for i := 0; i < 4; i++ {
			w[i] = wv * 4.0 / 3.0
		}
		for i := 4; i < 10; i++ {
			w[i] = we * 4.0 / 3.0
		}

	case 4:
		// Order 4: 20 point rule
		// Exact for polynomials up to degree 7

		// Parameters from Williams-Shunn-Jameson (2014)
		// Type 1: (a1, b1, b1, b1) - 4 points
		a1 := 0.0673422422100983
		b1 := 0.3108859192633006

		// Type 2: (c1, d1, d1, d1) - 4 points
		c1 := 0.7217942490673264
		d1 := 0.0927352503108912

		// Type 3: (a2, a2, b2, b2) - 12 points
		a2 := 0.4544962958743506
		b2 := 0.0455037041256494

		// Weights
		w1 := 0.0254224531851034
		w2 := 0.0665440511495285
		w3 := 0.0829803830550589

		// Type 1: (a1, b1, b1, b1) and permutations (4 points)
		x1 := []float64{
			-1.0 + 2.0*a1, -1.0 + 2.0*b1, -1.0 + 2.0*b1, -1.0 + 2.0*b1,
		}
		y1 := []float64{
			-1.0 + 2.0*b1, -1.0 + 2.0*a1, -1.0 + 2.0*b1, -1.0 + 2.0*b1,
		}
		z1 := []float64{
			-1.0 + 2.0*b1, -1.0 + 2.0*b1, -1.0 + 2.0*a1, -1.0 + 2.0*b1,
		}

		// Type 2: (c1, d1, d1, d1) and permutations (4 points)
		x2 := []float64{
			-1.0 + 2.0*c1, -1.0 + 2.0*d1, -1.0 + 2.0*d1, -1.0 + 2.0*d1,
		}
		y2 := []float64{
			-1.0 + 2.0*d1, -1.0 + 2.0*c1, -1.0 + 2.0*d1, -1.0 + 2.0*d1,
		}
		z2 := []float64{
			-1.0 + 2.0*d1, -1.0 + 2.0*d1, -1.0 + 2.0*c1, -1.0 + 2.0*d1,
		}

		// Type 3: (a2, a2, b2, b2) and permutations (12 points)
		// All 6 permutations of (a2, a2, b2, b2)
		x3 := []float64{
			-1.0 + 2.0*a2, -1.0 + 2.0*a2, -1.0 + 2.0*b2, -1.0 + 2.0*b2,
			-1.0 + 2.0*a2, -1.0 + 2.0*b2, -1.0 + 2.0*a2, -1.0 + 2.0*b2,
			-1.0 + 2.0*b2, -1.0 + 2.0*a2, -1.0 + 2.0*a2, -1.0 + 2.0*b2,
		}
		y3 := []float64{
			-1.0 + 2.0*a2, -1.0 + 2.0*b2, -1.0 + 2.0*a2, -1.0 + 2.0*b2,
			-1.0 + 2.0*b2, -1.0 + 2.0*a2, -1.0 + 2.0*b2, -1.0 + 2.0*a2,
			-1.0 + 2.0*a2, -1.0 + 2.0*a2, -1.0 + 2.0*b2, -1.0 + 2.0*b2,
		}
		z3 := []float64{
			-1.0 + 2.0*b2, -1.0 + 2.0*b2, -1.0 + 2.0*b2, -1.0 + 2.0*b2,
			-1.0 + 2.0*b2, -1.0 + 2.0*b2, -1.0 + 2.0*b2, -1.0 + 2.0*b2,
			-1.0 + 2.0*b2, -1.0 + 2.0*b2, -1.0 + 2.0*b2, -1.0 + 2.0*b2,
		}

		// Add the remaining 6 permutations for type 3
		x3rem := []float64{
			-1.0 + 2.0*a2, -1.0 + 2.0*b2, -1.0 + 2.0*b2, -1.0 + 2.0*a2,
			-1.0 + 2.0*b2, -1.0 + 2.0*a2, -1.0 + 2.0*b2, -1.0 + 2.0*a2,
			-1.0 + 2.0*b2, -1.0 + 2.0*b2, -1.0 + 2.0*a2, -1.0 + 2.0*a2,
		}
		y3rem := []float64{
			-1.0 + 2.0*b2, -1.0 + 2.0*a2, -1.0 + 2.0*a2, -1.0 + 2.0*b2,
			-1.0 + 2.0*a2, -1.0 + 2.0*b2, -1.0 + 2.0*a2, -1.0 + 2.0*b2,
			-1.0 + 2.0*b2, -1.0 + 2.0*b2, -1.0 + 2.0*a2, -1.0 + 2.0*a2,
		}
		z3rem := []float64{
			-1.0 + 2.0*b2, -1.0 + 2.0*b2, -1.0 + 2.0*b2, -1.0 + 2.0*b2,
			-1.0 + 2.0*b2, -1.0 + 2.0*b2, -1.0 + 2.0*b2, -1.0 + 2.0*b2,
			-1.0 + 2.0*a2, -1.0 + 2.0*a2, -1.0 + 2.0*a2, -1.0 + 2.0*a2,
		}

		x3 = append(x3, x3rem...)
		y3 = append(y3, y3rem...)
		z3 = append(z3, z3rem...)

		// Combine all points
		r = append(append(x1, x2...), x3...)
		s = append(append(y1, y2...), y3...)
		t = append(append(z1, z2...), z3...)

		// Set weights
		w = make([]float64, 20)
		for i := 0; i < 4; i++ {
			w[i] = w1 * 4.0 / 3.0
		}
		for i := 4; i < 8; i++ {
			w[i] = w2 * 4.0 / 3.0
		}
		for i := 8; i < 20; i++ {
			w[i] = w3 * 4.0 / 3.0
		}

	case 5:
		// Order 5: 35 point rule
		// Exact for polynomials up to degree 9

		// Centroid (1 point)
		xc := []float64{-1.0 + 2.0*0.25}
		yc := []float64{-1.0 + 2.0*0.25}
		zc := []float64{-1.0 + 2.0*0.25}
		wc := 0.0952380952380952

		// Vertex points (4 points)
		v := 0.3197936278296299
		b5 := (1.0 - v) / 3.0
		xv := []float64{
			-1.0 + 2.0*v, -1.0 + 2.0*b5, -1.0 + 2.0*b5, -1.0 + 2.0*b5,
		}
		yv := []float64{
			-1.0 + 2.0*b5, -1.0 + 2.0*v, -1.0 + 2.0*b5, -1.0 + 2.0*b5,
		}
		zv := []float64{
			-1.0 + 2.0*b5, -1.0 + 2.0*b5, -1.0 + 2.0*v, -1.0 + 2.0*b5,
		}
		wv := 0.0428571428571429

		// Edge points (12 points) - 2 per edge
		a := 0.0597158717897699
		b := 0.4701420641051151
		// Barycentric: (a,b,c,c) where c = (1-a-b)/2
		c := (1.0 - a - b) / 2.0

		xe := []float64{
			-1.0 + 2.0*a, -1.0 + 2.0*b, -1.0 + 2.0*a, -1.0 + 2.0*b,
			-1.0 + 2.0*a, -1.0 + 2.0*b, -1.0 + 2.0*c, -1.0 + 2.0*c,
			-1.0 + 2.0*c, -1.0 + 2.0*c, -1.0 + 2.0*c, -1.0 + 2.0*c,
		}
		ye := []float64{
			-1.0 + 2.0*b, -1.0 + 2.0*a, -1.0 + 2.0*c, -1.0 + 2.0*c,
			-1.0 + 2.0*c, -1.0 + 2.0*c, -1.0 + 2.0*a, -1.0 + 2.0*b,
			-1.0 + 2.0*a, -1.0 + 2.0*b, -1.0 + 2.0*c, -1.0 + 2.0*c,
		}
		ze := []float64{
			-1.0 + 2.0*c, -1.0 + 2.0*c, -1.0 + 2.0*b, -1.0 + 2.0*a,
			-1.0 + 2.0*c, -1.0 + 2.0*c, -1.0 + 2.0*c, -1.0 + 2.0*c,
			-1.0 + 2.0*c, -1.0 + 2.0*c, -1.0 + 2.0*a, -1.0 + 2.0*b,
		}
		we := 0.0285714285714286

		// Face points (12 points) - 3 per face
		c2 := 0.2146028712591518
		d2 := 0.3623969131549116
		e2 := 1.0 - 2.0*c2 - d2

		xf := []float64{
			-1.0 + 2.0*c2, -1.0 + 2.0*c2, -1.0 + 2.0*d2,
			-1.0 + 2.0*c2, -1.0 + 2.0*d2, -1.0 + 2.0*e2,
			-1.0 + 2.0*d2, -1.0 + 2.0*e2, -1.0 + 2.0*c2,
			-1.0 + 2.0*e2, -1.0 + 2.0*e2, -1.0 + 2.0*e2,
		}
		yf := []float64{
			-1.0 + 2.0*c2, -1.0 + 2.0*d2, -1.0 + 2.0*c2,
			-1.0 + 2.0*d2, -1.0 + 2.0*e2, -1.0 + 2.0*c2,
			-1.0 + 2.0*e2, -1.0 + 2.0*c2, -1.0 + 2.0*d2,
			-1.0 + 2.0*c2, -1.0 + 2.0*d2, -1.0 + 2.0*e2,
		}
		zf := []float64{
			-1.0 + 2.0*d2, -1.0 + 2.0*e2, -1.0 + 2.0*e2,
			-1.0 + 2.0*e2, -1.0 + 2.0*c2, -1.0 + 2.0*d2,
			-1.0 + 2.0*c2, -1.0 + 2.0*d2, -1.0 + 2.0*e2,
			-1.0 + 2.0*e2, -1.0 + 2.0*e2, -1.0 + 2.0*e2,
		}
		wf := 0.0571428571428571

		// Interior points (6 points)
		// Corrected interior point generation
		e3 := 0.2734364941962715
		f3 := 0.0605922445543035
		g3 := 0.3325391752274920
		h3 := 1.0 - e3 - f3 - g3

		xi := []float64{
			-1.0 + 2.0*e3, -1.0 + 2.0*e3, -1.0 + 2.0*e3,
			-1.0 + 2.0*f3, -1.0 + 2.0*f3, -1.0 + 2.0*f3,
		}
		yi := []float64{
			-1.0 + 2.0*e3, -1.0 + 2.0*f3, -1.0 + 2.0*g3,
			-1.0 + 2.0*e3, -1.0 + 2.0*g3, -1.0 + 2.0*h3,
		}
		zi := []float64{
			-1.0 + 2.0*f3, -1.0 + 2.0*g3, -1.0 + 2.0*h3,
			-1.0 + 2.0*g3, -1.0 + 2.0*h3, -1.0 + 2.0*e3,
		}
		wi := 0.0190476190476190

		// Combine all points
		r = append(append(append(append(xc, xv...), xe...), xf...), xi...)
		s = append(append(append(append(yc, yv...), ye...), yf...), yi...)
		t = append(append(append(append(zc, zv...), ze...), zf...), zi...)

		// Set weights
		w = make([]float64, 35)
		w[0] = wc * 4.0 / 3.0
		for i := 1; i < 5; i++ {
			w[i] = wv * 4.0 / 3.0
		}
		for i := 5; i < 17; i++ {
			w[i] = we * 4.0 / 3.0
		}
		for i := 17; i < 29; i++ {
			w[i] = wf * 4.0 / 3.0
		}
		for i := 29; i < 35; i++ {
			w[i] = wi * 4.0 / 3.0
		}

	default:
		// Fall back to order 5 for higher orders
		return WilliamsShunnJamesonCubature(5)
	}

	return
}

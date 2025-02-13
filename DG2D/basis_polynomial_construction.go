package DG2D

type ConstantVector struct {
	v [2]float64
}

func NewConstantVector(v1, v2 float64) *ConstantVector {
	return &ConstantVector{
		v: [2]float64{v1, v2},
	}
}

func (cv *ConstantVector) Eval() [2]float64 {
	return cv.v
}

func (cv *ConstantVector) Dot(f [2]float64) (dot float64) {
	dot = cv.v[0]*f[0] + cv.v[1]*f[1]
	return
}

func (cv *ConstantVector) Project(psi float64) (v [2]float64) {
	v = [2]float64{psi * cv.v[0], psi * cv.v[1]}
	return
}

type BaseVector struct {
	Eval       func(r, s float64) (v [2]float64)
	Dot        func(r, s float64, f [2]float64) (dot float64)
	Project    func(r, s float64, psi float64) (v [2]float64) // scalar mult psi
	Divergence func(r, s float64) (div float64)               // Div of vector
}

type PolynomialMultiplier struct {
	// This multiplies the VectorBase to produce a term
	Eval     func(r, s float64) (val float64)
	Gradient func(r, s float64) (grad [2]float64)
}
type VectorFunction struct {
	PolyMultiplier PolynomialMultiplier
	VectorBase     BaseVector
}

func (pt VectorFunction) Eval(r, s float64) (v [2]float64) {
	v = pt.VectorBase.Project(r, s, pt.PolyMultiplier.Eval(r, s))
	return
}
func (pt VectorFunction) Dot(r, s float64, b [2]float64) (dot float64) {
	v := pt.Eval(r, s)
	dot = b[0]*v[0] + b[1]*v[1]
	return
}

func (pt VectorFunction) Divergence(r, s float64) (div float64) {
	var (
		polyEval    = pt.PolyMultiplier.Eval(r, s)
		divBasis    = pt.VectorBase.Divergence(r, s)
		gradPoly    = pt.PolyMultiplier.Gradient(r, s)
		basisVector = pt.VectorBase
	)
	//    Div dot [P * e1, P * e2]
	//    =  (dP/dr)*e1 + P*(d(e1)/dr) + (dP/ds)*e2 + P*(d(e2)/ds)
	//    =  P * (d(e1)/dr + d(e2)/ds) + (dP/dr) * e1 + (dP/ds) * e2
	//    =  P * div([E]) + ([E])_dot_grad(P)
	div = polyEval*divBasis + basisVector.Dot(r, s, gradPoly)
	return
}

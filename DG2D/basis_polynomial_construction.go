package DG2D

type BasisVectorStruct struct {
	Eval       func(r, s float64) (v [2]float64)
	Dot        func(r, s float64, f [2]float64) (dot float64)
	Project    func(r, s float64, psi float64) (v [2]float64) // scalar mult psi
	Divergence func(r, s float64) (div float64)               // Div of vector
}

type BasisPolynomialMultiplier struct {
	// This multiplies the BasisVector to produce a term
	Eval     func(r, s float64) (val float64)
	Gradient func(r, s float64) (grad [2]float64)
}
type BasisPolynomialTerm struct {
	PolyMultiplier BasisPolynomialMultiplier
	BasisVector    BasisVectorStruct
}

func (pt BasisPolynomialTerm) Eval(r, s float64) (v [2]float64) {
	v = pt.BasisVector.Project(r, s, pt.PolyMultiplier.Eval(r, s))
	return
}
func (pt BasisPolynomialTerm) Dot(r, s float64, b [2]float64) (dot float64) {
	v := pt.Eval(r, s)
	dot = b[0]*v[0] + b[1]*v[1]
	return
}

func (pt BasisPolynomialTerm) Divergence(r, s float64) (div float64) {
	var (
		polyEval    = pt.PolyMultiplier.Eval(r, s)
		divBasis    = pt.BasisVector.Divergence(r, s)
		gradPoly    = pt.PolyMultiplier.Gradient(r, s)
		basisVector = pt.BasisVector
	)
	//    Div dot [P * e1, P * e2]
	//    =  (dP/dr)*e1 + P*(d(e1)/dr) + (dP/ds)*e2 + P*(d(e2)/ds)
	//    =  P * (d(e1)/dr + d(e2)/ds) + (dP/dr) * e1 + (dP/ds) * e2
	//    =  P * div([E]) + ([E])_dot_grad(P)
	div = polyEval*divBasis + basisVector.Dot(r, s, gradPoly)
	return
}

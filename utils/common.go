package utils

const (
	NODETOL = 1.e-12
)

type EvalOp uint8

const (
	Equal EvalOp = iota
	Less
	Greater
	LessOrEqual
	GreaterOrEqual
)

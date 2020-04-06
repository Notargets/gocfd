package utils

type EvalOp uint8

const (
	Equal EvalOp = iota
	Less
	Greater
	LessOrEqual
	GreaterOrEqual
)

type MathOp uint8

const (
	Abs MathOp = iota
	Sqrt
	Pow
)

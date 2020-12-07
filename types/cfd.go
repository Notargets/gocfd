package types

//go:generate stringer -type=BCFLAG

type BCFLAG uint8

const (
	BC_None BCFLAG = iota
	BC_In
	BC_Dirichlet
	BC_Slip
	BC_Far
	BC_Wall
	BC_Cyl
	BC_Neuman
	BC_Out
	BC_IVortex
)

var BCNameMap = map[string]BCFLAG{
	"inflow":    BC_In,
	"in":        BC_In,
	"out":       BC_Out,
	"outflow":   BC_Out,
	"wall":      BC_Wall,
	"far":       BC_Far,
	"cyl":       BC_Cyl,
	"dirichlet": BC_Dirichlet,
	"neuman":    BC_Neuman,
	"slip":      BC_Slip,
}

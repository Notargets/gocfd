package utils

import "strings"

// BCType represents boundary condition types for CFD/MHD simulations
type BCType uint16

// Boundary condition constants for typical CFD/MHD applications
const (
	// BCNone indicates no boundary condition (interior face)
	BCNone BCType = iota

	// Flow boundary conditions
	BCInflow   // Inflow/inlet boundary
	BCOutflow  // Outflow/outlet boundary
	BCWall     // No-slip wall
	BCSlipWall // Slip/inviscid wall
	BCSymmetry // Symmetry plane
	BCPeriodic // Periodic boundary
	BCFarfield // Far-field boundary

	// Thermal boundary conditions
	BCIsothermal // Fixed temperature
	BCAdiabatic  // No heat flux
	BCHeatFlux   // Prescribed heat flux

	// Electromagnetic boundary conditions (for MHD)
	BCPerfectConductor // Perfect electric conductor
	BCInsulator        // Perfect insulator
	BCMagneticWall     // Perfectly conducting magnetic wall

	// Mathematical boundary conditions
	BCDirichlet // Fixed value
	BCNeumann   // Fixed gradient/flux
	BCRobin     // Mixed/Robin condition

	// Special boundary conditions
	BCInterface      // Interface between domains
	BCMovingWall     // Moving wall boundary
	BCPressureOutlet // Pressure-specified outlet
	BCMassFlowInlet  // Mass flow rate inlet
	BCVelocityInlet  // Velocity-specified inlet

	// Parallel/domain decomposition
	BCPartitionBoundary // Boundary between parallel partitions

	// User-defined (reserve space for custom BCs)
	BCUserDefined1
	BCUserDefined2
	BCUserDefined3
	BCUserDefined4
	BCUserDefined5
)

// String returns the string representation of a BCType
func (bc BCType) String() string {
	names := map[BCType]string{
		BCNone:              "None",
		BCInflow:            "Inflow",
		BCOutflow:           "Outflow",
		BCWall:              "Wall",
		BCSlipWall:          "SlipWall",
		BCSymmetry:          "Symmetry",
		BCPeriodic:          "Periodic",
		BCFarfield:          "Farfield",
		BCIsothermal:        "Isothermal",
		BCAdiabatic:         "Adiabatic",
		BCHeatFlux:          "HeatFlux",
		BCPerfectConductor:  "PerfectConductor",
		BCInsulator:         "Insulator",
		BCMagneticWall:      "MagneticWall",
		BCDirichlet:         "Dirichlet",
		BCNeumann:           "Neumann",
		BCRobin:             "Robin",
		BCInterface:         "Interface",
		BCMovingWall:        "MovingWall",
		BCPressureOutlet:    "PressureOutlet",
		BCMassFlowInlet:     "MassFlowInlet",
		BCVelocityInlet:     "VelocityInlet",
		BCPartitionBoundary: "PartitionBoundary",
		BCUserDefined1:      "UserDefined1",
		BCUserDefined2:      "UserDefined2",
		BCUserDefined3:      "UserDefined3",
		BCUserDefined4:      "UserDefined4",
		BCUserDefined5:      "UserDefined5",
	}

	if name, ok := names[bc]; ok {
		return name
	}
	return "Unknown"
}

// BCNameMap provides a mapping from common boundary condition names to BCType
// Keys are lowercase for case-insensitive matching
// This can be extended by applications to support mesh-specific naming conventions
var BCNameMap = map[string]BCType{
	// Common variations for inflow
	"inlet":           BCInflow,
	"inflow":          BCInflow,
	"velocity_inlet":  BCVelocityInlet,
	"mass_flow_inlet": BCMassFlowInlet,

	// Common variations for outflow
	"outlet":          BCOutflow,
	"outflow":         BCOutflow,
	"exit":            BCOutflow,
	"pressure_outlet": BCPressureOutlet,

	// Wall variations
	"wall":          BCWall,
	"no_slip":       BCWall,
	"noslip":        BCWall,
	"slip":          BCSlipWall,
	"slip_wall":     BCSlipWall,
	"inviscid_wall": BCSlipWall,

	// Other boundaries
	"symmetry":   BCSymmetry,
	"symmetric":  BCSymmetry,
	"farfield":   BCFarfield,
	"far_field":  BCFarfield,
	"freestream": BCFarfield,
	"periodic":   BCPeriodic,

	// Thermal
	"isothermal": BCIsothermal,
	"adiabatic":  BCAdiabatic,
	"heat_flux":  BCHeatFlux,

	// Mathematical
	"dirichlet": BCDirichlet,
	"neumann":   BCNeumann,
	"robin":     BCRobin,

	// Interface
	"interface": BCInterface,
	"internal":  BCInterface,
}

// To add custom BC names, applications can modify BCNameMap:
//   utils.BCNameMap["my_custom_bc"] = utils.BCUserDefined1

// ParseBCName converts a boundary condition name string to BCType
// The matching is case-insensitive and trims whitespace
func ParseBCName(name string) BCType {
	// Convert to lowercase for case-insensitive matching
	lowerName := strings.ToLower(strings.TrimSpace(name))

	if bcType, ok := BCNameMap[lowerName]; ok {
		return bcType
	}
	// Default to wall for unknown types
	return BCWall
}


# Quick Start Guide


## Building on Ubuntu Linux
Ensure you have the Go language installed and available in your PATH. Then, install the necessary prerequisites:
```bash

sudo apt update
sudo apt install libx11-dev libxi-dev libxcursor-dev libxrandr-dev libxinerama-dev mesa-common-dev libgl1-mesa-dev libxxf86vm-dev
make
```

## Running Test Cases
```

### 1D Shock Tube Test Case
#### Without Graphics
gocfd 1D

#### With Graphics
export DISPLAY=:0
gocfd 1D -g

### 2D Airfoil Test Case
gocfd 2D -F test_cases/Euler2D/naca_12/mesh/nacaAirfoil-base.su2 -I test_cases/Euler2D/naca_12/input-wall.yaml -g -z 0.08
```

## Code Review Guide

To review the physics implementation, explore the code within the `model_problems` directory. Each file implements a physical model or an additional numerical method.

For insights into the math/matrix library, refer to `utils/matrix_extended.go` and `utils/vector_extended.go`. A chained operator syntax is used here to prioritize reuse, reduce copying, and clarify operand dimensionality. Although not exhaustive, especially regarding value assignment and indexing, it mirrors and implements the functionalities found in related texts and MATLAB, thus proving functional and useful. The syntax resembles reverse Polish notation (RPN), where values accumulate through chained operations.

For example, consider the following implementation:
```rhsE = - Dr * FluxH .* Rx ./ epsilon```. Here, `Dr` stands for the "Derivative Matrix", `Rx` represents (1 / J) applying the transform of R to X, and `epsilon` is the metallic impedance, all applied to the Flux matrix:
```
	RHSE = el.Dr.Mul(FluxH).ElMul(el.Rx).ElDiv(c.Epsilon).Scale(-1)
```

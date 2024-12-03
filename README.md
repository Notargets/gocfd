# gocfd
A CFD Solver implemented in the Go programming language.

Implements the Direct Flux Reconstruction (DFR) method within Galerkin Discontinuous Finite Elements (GDFE) for unstructured meshes.

For detailed progress logs, see [NOTES_Index.md.](NOTES_Index.md). Currently, there is a functional solver for 2D Euler Equations allowing the input of unstructured meshes and boundary conditions to simulate flowfields. Effectively capturing shock waves and contact discontinuities is underway.

The next steps, post shock capture capability, will involve incorporating viscous flowfields and likely adding various implicit time integration solvers.

| NACA 0012 Airfoil at M=0.3, Alpha=6, Roe flux, Local Time Stepping | M=0.5, Alpha=0, Roe Flux, 1482 O(2) Elements, Converged |
|:------------------------------------------------------------------:|--------------------------------------------------------:|
|               ![](images/naca12_2d_m0.3_a6_roe.gif)                |                 ![](images/naca12_2d_m0.5_aoa0_Roe.PNG) |

|                           Density                            |                            X Momentum                             |                  Density                   |
|:------------------------------------------------------------:|:-----------------------------------------------------------------:|:------------------------------------------:|
| ![](images/render-mesh-isentropic-vortex-initial-zoom-7.PNG) | ![](images/render-mesh-isentropic-vortex-initial-zoom-7-rhoU.png) | ![](images/vortex-1-2-4-7-lax-cropped.gif) |


## Discontinuous Galerkin Method for solving systems of equations - CFD, CEM, ... hydrodynamics-fusion (simulate the Sun), etc!

##### Credits to:
- Jan S. Hesthaven and Tim Warburton for their excellent text "Nodal Discontinuous Galerkin Methods" (2007)
- J. Romero, K. Asthana, and Antony Jameson for "A Simplified Formulation of the Flux Reconstruction Method" (2015) for the
  DFR approach with Raviart-Thomas elements

### Objectives

1) Implement a complete 3D solver for unstructured CFD (and possibly MHD) using the Discontinuous Galerkin (DG) method
2) Optimize for GPUs and groups of GPUs, taking advantage of the nature of the natural parallelism of DG methods
3) Prove the accuracy of the CFD solver for predicting flows with turbulence, shear flows and strong temperature gradients
4) Make the solver available for use as an open source tool

It is important to me that the code implementing the solver be as simple as possible so that it can be further developed
and extended. There are other projects that have achieved some of the above, most notably the
HiFiLES project (formerly at https://hifiles.stanford.edu/) project, which has demonstrated high accuracy for turbulence 
problems and some transonic flows with shock waves and is open source. I personally find that C++ code is very difficult 
to understand due to the heavy usage of indirection and abstraction, which makes an already complex subject unnecessarily
more difficult. I feel that the Go language makes it easier to develop straightforward, more easily understandable code
paths, while providing similar if not equivalent optimality and higher development efficiency than C++.

### Why do this work?

I studied CFD in graduate school in 1987 and worked for Northrop for 10 years building and using CFD methods to design
and debug airplanes and propulsion systems. During my time applying CFD, I had some great success and some notable
failures in getting useful results from the CFD analysis. The most common theme in the failures: flows with thermal
gradients, shear flows, and vortices were handled very poorly by all known usable Finite Volume methods.

Then, last year (2019), I noticed there were some amazing looking results appearing on YouTube and elsewhere showing
well-resolved turbulent eddies and shear flows using this new "Discontinuous Galerkin Finite Elements" method...

====
# Quick Start Guide

## Building on Ubuntu Linux

First, ensure the Go language is installed and available in your PATH. Proceed to install the necessary prerequisites:

```
sudo apt update
sudo apt install libx11-dev libxi-dev libxcursor-dev libxrandr-dev libxinerama-dev mesa-common-dev libgl1-mesa-dev libxxf86vm-dev

make
```
```

## Running Test Cases

### 1D Shock Tube Test Case
#### Without Graphics
gocfd 1D

#### With Graphics Enabled
export DISPLAY=:0
gocfd 1D -g

### 2D Airfoil Test Case
gocfd 2D -F test_cases/Euler2D/naca_12/mesh/nacaAirfoil-base.su2 -I test_cases/Euler2D/naca_12/input-wall.yaml -g -s 50 -z 0.08
```

## Code Review Guidelines

To review the physics implementation, explore the code within the `model_problems` directory. Each file implements a physical model or an additional numerical method.

For insights into the math/matrix library, see `utils/matrix_extended.go` and `utils/vector_extended.go`. They utilize a chained operator syntax prioritizing reuse, reducing copying, and clarifying operand dimensionality. While not exhaustive, especially in terms of value assignment and indexing, this approach mirrors the functionalities found in related texts and MATLAB, proving both functional and useful. The syntax resembles reverse Polish notation (RPN), where values accumulate through chained operations.

For example, consider this implementation:
```rhsE = - Dr * FluxH .* Rx ./ epsilon```. Here, `Dr` stands for the "Derivative Matrix", `Rx` represents (1 / J) applying the transform of R to X, and `epsilon` is the metallic impedance, all applied to the Flux matrix:
```
RHSE = el.Dr.Mul(FluxH).ElMul(el.Rx).ElDiv(c.Epsilon).Scale(-1)
```

## Requirements to run the code
Here is what I'm using as a platform:
```
me@home:bash# go version
go version go1.23.3 linux/amd64

me@home:bash# lsb_release -a
No LSB modules are available.
Distributor ID:	Ubuntu
Description:	Ubuntu 24.04.1 LTS
Release:	24.04
Codename:	noble
...
```
You also need to install some X11 and OpenGL related packages in Ubuntu, like this:
```
apt update
apt install libx11-dev libxi-dev libxcursor-dev libxrandr-dev libxinerama-dev mesa-common-dev libgl1-mesa-dev
```
A proper build should go like this:
```
me@home:bash# make
go fmt ./...  && go install ./...
go: downloading gonum.org/v1/gonum v0.7.0
go: downloading github.com/notargets/avs v0.0.7-0.20201217183319-0f9c8f3d02f3
go: downloading github.com/james-bowman/sparse v0.0.0-20200204164517-b588421ac5da
go: downloading github.com/spf13/cobra v1.0.0
go: downloading github.com/spf13/viper v1.7.0
go: downloading github.com/mitchellh/go-homedir v1.1.0
go: downloading github.com/ghodss/yaml v1.0.0
go: downloading github.com/pkg/profile v1.5.0
go: downloading gopkg.in/yaml.v2 v2.3.0
go: downloading github.com/hashicorp/hcl v1.0.0
go: downloading github.com/pelletier/go-toml v1.8.0
go: downloading github.com/spf13/afero v1.2.2
go: downloading github.com/fsnotify/fsnotify v1.4.9
go: downloading github.com/magiconair/properties v1.8.1
go: downloading github.com/mitchellh/mapstructure v1.3.2
go: downloading github.com/spf13/cast v1.3.1
go: downloading github.com/spf13/jwalterweatherman v1.1.0
go: downloading github.com/spf13/pflag v1.0.5
go: downloading github.com/subosito/gotenv v1.2.0
go: downloading gopkg.in/ini.v1 v1.57.0
go: downloading golang.org/x/sys v0.0.0-20220722155257-8c9f86f7a55f
go: downloading golang.org/x/text v0.3.8
go: downloading github.com/go-gl/gl v0.0.0-20190320180904-bf2b1f2f34d7
go: downloading github.com/go-gl/glfw/v3.3/glfw v0.0.0-20200420212212-258d9bec320e
go: downloading github.com/gonum/floats v0.0.0-20181209220543-c233463c7e82
go: downloading github.com/gonum/internal v0.0.0-20181124074243-f884aa714029
go: downloading github.com/go-gl/glfw v0.0.0-20190409004039-e6da0acd62b1
run this -> $GOPATH/bin/gocfd
```

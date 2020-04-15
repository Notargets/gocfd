# gocfd
Awesome CFD solver written in Go 

## An implementation of the Discontinuous Galerkin Method for solving systems of equations

## Requirements
Here is what I'm using:
```
me@home:bash# go version
go version go1.14 linux/amd64
me@home:bash# cat /etc/os-release
NAME="Ubuntu"
VERSION="18.04.4 LTS (Bionic Beaver)"
...
```
A proper build should go like this:
```
me@home:bash# $ make
go fmt ./...  && go install ./...
run this -> $GOPATH/bin/gocfd
me@home:bash# /gocfd$ gocfd --help
Usage of gocfd:
  -K int
        Number of elements in model (default 60)
  -N int
        polynomial degree (default 8)
  -delay int
        milliseconds of delay for plotting (default 0)
  -graph
        display a graph while computing solution
  -model int
        model to run: 0 = Advect1D, 1 = Maxwell1D (default 1)
```

### Model Problem Example #2: Maxwell's Equations solved in a 1D Cavity

The Maxwell equations are solved in a 1D metal cavity with a change of material half way through the domain. The initial condition is a sine wave for the E (electric) field in the left half of the domain, and zero for E and H everywhere else. The E field is zero on the boundary (face flux out = face flux in) and the H field passes through unchanged (face flux zero), corresponding to a metallic boundary.



Run the example with graphics like this:
```
bash# make
bash# gocfd -model 1 -delay 0 -graph -K 80 -N 5
```

Unlike the advection equation model problem, this solver does have unstable points in the space of K (element count) and N (polynomial degree). So far, it appears that the polynomial degree must be >= 5 for stability, otherwise aliasing occurs, where even/odd modes are excited among grid points.

In the example pictured, there are 80 elements (K=80) and the element polynomial degree is 5 (N=5).

#### Initial State
![](images/Maxwell1D-cavity0.PNG)

#### Intermediate State
![](images/Maxwell1D-cavity.PNG)

#### First Mode
![](images/Maxwell1D-cavity2.PNG)

#### Second Mode
![](images/Maxwell1D-cavity4.PNG)

### Model Problem Example #1: Advection Equation
<span style="display:block;text-align:center">![](images/Advect1D-0.PNG)</span>

The first model problem is 1D Advection with a left boundary driven sine wave. You can run it with integrated graphics like this:
```
bash# gocfd -model 0 -delay 0 -graph -K 80 -N 5
```

In the example pictured, there are 80 elements (K=80) and the element polynomial degree is 5 (N=5).

<span style="display:block;text-align:center">![](images/Advect1D-1.PNG)</span>

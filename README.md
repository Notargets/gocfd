# gocfd
Awesome CFD solver written in Go 

## An implementation of the Discontinuous Galerkin Method for solving systems of equations

### Credits to Jan S. Hesthaven and Tim Warburton for their excellent text "Nodal Discontinuous Galerkin Methods" (2007)

### QuickStart

Using Ubuntu Linux, do the following first:
```
me@home:bash# sudo apt update
me@home:bash# sudo apt install libx11-dev libxi-dev libxcursor-dev libxrandr-dev libxinerama-dev mesa-common-dev libgl1-mesa-dev
me@home:bash# make
with graphics:
me@home:bash# export DISPLAY=:0
me@home:bash# gocfd -graph
without graphics:
me@home:bash# gocfd
```

### Current Status (May 9, 2020): Direct Flux Reconstruction implemented for 1D Advection, moving to implement for the other two 1D model problems

DFR works for Advection (-model 3) and seems to improve accuracy and physicality.

The primary differences between using DFR and traditional Nodal Discontinuous Galerkin:
1) Instead of using only the primitive variables in the right hand side, we use the Flux directly for a hyperbolic problem
2) The flux is a globally composed variable and is differentiated after being "reconstructed" to be C(1)
3) The reconstruction technique follows Jameson(2014) in that we coerce the flux values of each face to be consistent values at each face, then use the same Np Gauss-Lobato nodes and metrics to develop the derivative(s) of flux

For (3) in this case, I used the same Lax-Friedrichs flux calculation from the text, then averaged the values from each side of each face together to make a single consistent face value shared by neighboring elements.

### (May 1, 2020): Researching Direct Flux Reconstruction methods
During testing of the method in 1D as outlined in the text, it became clear that the slope limiter is quite crude and is degrading the physicality of the solution. The authors were clear that this is just for example use, and now I'm convinced of it!

My first round of research yielded the [Direct Flux Reconstruction](research/filters_and_flux_limiters/Romero2015.pdf) technique from a 2014 paper by Antony Jameson, et al (one of my favorite CFD people of all time). The technique is extremely simple and has the great promise of extending the degree of accuracy to the Flux terms of more complex equations, in addition to enabling the use of flux limiters that have been proven for flows with discontinuities.

In my past CFD experience, discontinuities in solving the Navier Stokes equations are not limited to shock waves. Rather, we find shear flows have discontinuities that are similar enough to shock waves that shock finding techniques used to guide the application of numerical diffusion can be active at the edge of boundary layers and at other critical regions, which often leads to inaccuracies far from shocks. The text's gradient based limiter would clearly suffer from this kind of problem among others.

## Requirements to run the code
Here is what I'm using as a platform:
```
me@home:bash# go version
go version go1.14 linux/amd64
me@home:bash# cat /etc/os-release
NAME="Ubuntu"
VERSION="18.04.4 LTS (Bionic Beaver)"
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
### Current Status: Researching Direct Flux Reconstruction methods

During testing of the method in 1D as outlined in the text, it became clear that the slope limiter is quite crude and is degrading the physicality of the solution. The authors were clear that this is just for example use, and now I'm convinced of it!

My first round of research yielded the [Direct Flux Reconstruction](research/filters_and_flux_limiters/Romero2015.pdf) technique from a 2014 paper by Antony Jameson, et al (one of my favorite CFD people of all time). The technique is extremely simple and has the great promise of extending the degree of accuracy to the Flux terms of more complex equations, in addition to enabling the use of flux limiters that have been proven for flows with discontinuities.

In my past CFD experience, discontinuities in solving the Navier Stokes equations are not limited to shock waves. Rather, we find shear flows have discontinuities that are similar enough to shock waves that shock finding techniques used to guide the application of numerical diffusion can be active at the edge of boundary layers and at other critical regions, which often leads to inaccuracies far from shocks. The text's gradient based limiter would clearly suffer from this kind of problem among others.

### Model Problem Example #3a: Euler's Equations in 1D - Shock Collision

This is an interesting problem because of the temperature remainder after the collision. In the plot, temperature is red, density is blue, and velocity is orange. After the shocks pass out of the domain, the remaining temperature "bubble" can't dissipate, because the Euler equations have no mechanism for temperature diffusion.

This case is obtained by initializing the tube as with the Sod tube, but leaving the exit boundary at the left side values (In == Out). This produces a left running shock wave that meets with the shock moving right.

#### T = 0.36, 1000 Elements
![](images/Euler1D-MidTube-K1000-N6-T.36.PNG)

### Model Problem Example #3: Euler's Equations in 1D - Sod's shock tube

The 1D Euler equations are solved with boundary and initial conditions for the Sod shock tube problem. There is an analytic solution for this case and it is widely used to test shock capturing ability of a solver.

Run the example with graphics like this:
```
bash# make
bash# gocfd -model 2 -graph -K 250 -N 1
```

You can also target a final time for the simulation using the "-FinalTime" flag. You will have to use CTRL-C to exit the simulation when it arrives at the target time. This leaves the plot on screen so you can screen cap it.
```
bash# gocfd -model 2 -graph -K 250 -N 1 -FinalTime 0.2
```
#### T = 0.2, 60 Elements
Linear Elements | 10th Order Elements
:-------------------------:|:-------------------------:
![](images/Euler1D-SOD-K60-N1-T0.2.PNG) | ![](images/Euler1D-SOD-K60-N10-T0.2.PNG)

#### T = 0.2, 250 Elements
Linear Elements | 10th Order Elements
:-------------------------:|:-------------------------:
![](images/Euler1D-SOD-K250-N1-T0.2.PNG) | ![](images/Euler1D-SOD-K250-N10-T0.2.PNG)

#### T = 0.2, 500 Elements
Linear Elements | 10th Order Elements
:-------------------------:|:-------------------------:
![](images/Euler1D-SOD-K500-N1-T0.2.PNG) | ![](images/Euler1D-SOD-K500-N10-T0.2.PNG)


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

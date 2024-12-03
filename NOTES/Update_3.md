## Update: Model Problem Example #3a: Euler's Equations in 1D - Shock Collision

This is an interesting problem because of the temperature remainder after the collision. In the plot, temperature is red, density is blue, and velocity is orange. After the shocks pass out of the domain, the remaining temperature "bubble" can't dissipate, because the Euler equations have no mechanism for temperature diffusion.

This case is obtained by initializing the tube as with the Sod tube, but leaving the exit boundary at the left side values (In == Out). This produces a left running shock wave that meets with the shock moving right.

#### T = 0.36, 1000 Elements
![](../images/Euler1D-MidTube-K1000-N6-T.36.PNG)

### Update: Model Problem Example #3: Euler's Equations in 1D - Sod's shock tube

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
|             Linear Elements             |           10th Order Elements            |
|:---------------------------------------:|:----------------------------------------:|
| ![](../images/Euler1D-SOD-K60-N1-T0.2.PNG) | ![](../images/Euler1D-SOD-K60-N10-T0.2.PNG) |

#### T = 0.2, 250 Elements
|             Linear Elements              |            10th Order Elements            |
|:----------------------------------------:|:-----------------------------------------:|
| ![](../images/Euler1D-SOD-K250-N1-T0.2.PNG) | ![](../images/Euler1D-SOD-K250-N10-T0.2.PNG) |

#### T = 0.2, 500 Elements
|             Linear Elements              |            10th Order Elements            |
|:----------------------------------------:|:-----------------------------------------:|
| ![](../images/Euler1D-SOD-K500-N1-T0.2.PNG) | ![](../images/Euler1D-SOD-K500-N10-T0.2.PNG) |


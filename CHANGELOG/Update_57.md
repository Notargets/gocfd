# Shock Capturing - Success with Polynomial Limiter!

## The problem was always the solution polynomial, so limit that

When we look at the problem of edge interpolation experiencing Gibbs 
oscillations, the core issue is really that the polynomial at the interior 
that represents the solution has too much energy overall. The edge 
interpolations are just an expression of that, so in the end we have to do 
something to tamp down the energy of the solution polynomial itself when it 
is unable to express the discontinuity. The artificial dissipation approach 
does this incrementally, but in the end I found it was insufficient to 
mitigate enough energy in the polynomial to remove the edge interpolation 
overshoots. When I tried increasing the dissipation too high, I ended up 
with either washed out solutions or instability. How about just modulating 
the solution polynomial when it is "tipped over"?

I had previously tried using the (now ancient) Barth Jesperson limiter, but 
had found that I couldn't remove the buzzing and so couldn't converge 
solutions with it. If we could get something that smoothly transitions to 
damping when needed, we should be able to get converged solutions without 
harming the solution accuracy away from the "troubled cells" where we have 
discontinuities. Ideally, we can have an accurate first order solution at 
the most shocked cells and not have nearby cells affected.

Below is what I found to be effective when used to limit the solution 
polynomial *at the last Runge Kutta stage*, before the solution is 
interpolated to the edges at that stage.

        Q_Limited = (1. - alpha) * Q + alpha * Q_Mean
        alpha = 1. - exp(-Beta * sigma)
            sigma is the shock indicator
            Beta is a constant, around 3.->10.

## Experiments with existing shock indicator and Persson viscocity

| Re-entry vehicle? Mach 5, AOA 35 degrees, P=2 |
|-----------------------------------------------|
| ![](naca0012-aoa35-M5.png)                    |

I can now run ridiculous configurations like the NACA0012 airfoil 
re-entering the atmosphere at 35 degrees angle of attack and Mach 5, at 
polynomial orders up to P=6. In experiments with the shock indicator, I've 
learned there are some structural issues and a bug in the current 
implementation that I'd like to refine, but this is definitely on the right 
track. It only took four years LOL.

I think the primary issue I've had in confronting this problem is that I've 
been unwilling to limit the solution itself, thinking that it's a betrayal 
of the high order accuracy approach. I've come to realize that when the 
polynomial itself is unable to capture the discontinuity, a model can be 
substituted for just the affected cells that still has physical fidelity. 
Using the current filter, the model is a hybrid of the calculated cell mean 
and the polynomial that is clearly able to suppress the polynomial, but it 
appears that the accuracy of the solution is accurate overall.

## Next steps

I'm going to refine the shock finder in three ways:
1) Fix a bug in the calculation of the trigger energy - it's currently P^-4 
   when it should be (P+1)^-4.
2) Change the Sin() ramp of the trigger to Persson's original power ratio
3) Fix a bug in the moment used to determine a cell is "in trouble" to more 
   accurately represent the energy and make it run faster
syms x y beta gamma t x0 y0 bt
r2=(x-t-x0)^2+(y-y0)^2
bt=beta*exp(1-r2)
bt2=bt^2
%isentropic vortex solution for Euler equations in 2D
u=1-bt*(y-y0)
v=bt*(x-x0)
rho=(1-bt2*((gamma-1)/gamma)/4)^(1/(gamma-1))
q=0.5*rho*(u^2+v^2)
p=rho^gamma
E=(p/(gamma-1)+q)/rho % generic calculation of E, not specific to this case, must happen after rho is defined
F = [ rho*u, rho*u^2+p, rho*v, u*(E+p) ]
G = [ rho*v, rho*u*v, rho*v^2+p, v*(E+p) ]
div = diff(F,x)+diff(G,y)

syms x y gamma
%2D Polynomial field
rho=x^2+y^2;
rhou=x^2;
rhov=y^2;
E=10*(x^2+y^2);
p = rho^gamma
u = rhou/rho;
v = rhov/rho;
U = [ rho, rhou, rhov, E];
F = [ rhou, rho*u^2+p, rho*u*v, u*(E+p) ];
G = [ rhov, rho*u*v, rho*v^2+p, v*(E+p) ];
div = diff(F,x)+diff(G,y);
fprintf('Code for Divergence of F and G Fluxes\n%s\n',ccode(div));
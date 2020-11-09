syms a b c d x y gamma
%2D Polynomial field
rho=a*abs(x)+b*abs(y);
u = c*x; v = d*y;
rhou=rho*u; rhov=rho*v;
p=rho^gamma;
q=0.5*rho*(u^2+v^2);
E=p/(gamma-1)+q;
U = [ rho, rhou, rhov, E];
F = [ rhou, rho*u^2+p, rho*u*v, u*(E+p) ];
G = [ rhov, rho*u*v, rho*v^2+p, v*(E+p) ];
div = diff(F,x)+diff(G,y);
fprintf('Code for Divergence of F and G Fluxes\n%s\n',ccode(div));
fprintf('Code for U \n%s\n%s\n%s\n%s\n',ccode(U));
plotRho(x,y) = subs(rho,[a,b],[1.,1.]);
plotDiv(x,y) = subs(div,[a,b,c,d,gamma],[1.,1.,1.,1.,1.4]);
plotU(x,y) = subs(U,[a,b,c,d,gamma],[1.,1.,1.,1.,1.4]);
subs(U(4),[a,b,c,d,gamma],[1.,1.,1.,1.,1.4])
x = linspace(-10,10,20);
y = linspace(-10,10,20);
[X,Y] = meshgrid(x,y);
Z = plotU(X,Y);
figure('Name','rho','NumberTitle','off');
contour(X,Y,double(Z{1}),10, 'ShowText','on');
figure('Name','rhoU','NumberTitle','off');
contour(X,Y,double(Z{2}),10, 'ShowText','on');
figure('Name','rhoV','NumberTitle','off');
contour(X,Y,double(Z{3}),10, 'ShowText','on');
figure('Name','E','NumberTitle','off');
contour(X,Y,double(Z{4}),10, 'ShowText','on');
divZ = plotDiv(X,Y);
figure('Name','Divergence of Rho','NumberTitle','off');
contour(X,Y,double(divZ{1}),10, 'ShowText','on');
figure('Name','Divergence of RhoU','NumberTitle','off');
contour(X,Y,double(divZ{2}),10, 'ShowText','on');
figure('Name','Divergence of RhoV','NumberTitle','off');
contour(X,Y,double(divZ{3}),10, 'ShowText','on');
figure('Name','Divergence of E','NumberTitle','off');
contour(X,Y,double(divZ{4}),10, 'ShowText','on');
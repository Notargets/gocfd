syms x y beta gamma t x0 y0
r2=(x-t-x0)^2+(y-y0)^2;
bt=beta*exp(1-r2)/(2*pi);
bt2=bt^2;
%isentropic vortex solution for Euler equations in 2D
u=1-bt*(y-y0);
v=bt*(x-x0);
rho=(1-bt2*(gamma-1)/(4*gamma))^(1/(gamma-1));
q=0.5*rho*(u^2+v^2);
p=rho^gamma;
E=(p/(gamma-1)+q)/rho; % generic calculation of E, not specific to this case, must happen after rho is defined
U = [ rho, rho*u, rho*v, rho*E];
F = [ rho*u, rho*u^2+p, rho*u*v, u*(E+p) ];
G = [ rho*v, rho*u*v, rho*v^2+p, v*(E+p) ];
div = diff(F,x)+diff(G,y);
fprintf('Code for Divergence of F and G Fluxes\n%s\n',ccode(div));
disp(ccode(div(1)));
disp(ccode(div(2)));
disp(ccode(div(3)));
disp(ccode(div(4)));

%Test for correctness: use a finite difference stencil approximation to the
%original function
Beta = 5; X0=5; Y0=0; T=0; GAMMA=1.4;
X=5; Y=0;
DX=0.001; DY=0.001;
args=[x,y,beta,gamma,t,x0,y0];
vals=[X-DX/2,Y,Beta,GAMMA,T,X0,Y0];
Fxm = subs(F,args,vals);
vals=[X+DX/2,Y,Beta,GAMMA,T,X0,Y0];
Fxp = subs(F,args,vals);
Fx=(Fxp-Fxm)/DX;

vals=[X,Y-DY/2,Beta,GAMMA,T,X0,Y0];
Gym = subs(G,args,vals);
vals=[X,Y+DY/2,Beta,GAMMA,T,X0,Y0];
Gyp = subs(G,args,vals);
Gy=(Gyp-Gym)/DY;
Div=Fx+Gy;
fprintf('Divergence at X=[%f,%f] = [%f,%f,%f,%f]\n',double(X), double(Y),double(Div));

vals=[X,Y,Beta,GAMMA,T,X0,Y0];
adiv=subs(div,args,vals);
fprintf('Analytical Divergence at X=[%f,%f] = [%f,%f,%f,%f]\n',double(X), double(Y), double(adiv));

residual = double(adiv-Div);
fprintf('*** SHOULD BE ZERO *** Residual at X=[%f,%f] = [%f,%f,%f,%f]\n',double(X), double(Y), double(residual));

fprintf('F at X=[%f,%f] = [%f,%f,%f,%f]\n',double(X), double(Y), double(subs(F,args,vals)));
fprintf('G at X=[%f,%f] = [%f,%f,%f,%f]\n',double(X), double(Y), double(subs(G,args,vals)));
fprintf('Code for Flux in F and G directions\n%s\n%s\n',ccode(F),ccode(G));

fprintf('rho,rhoU,rhoV,rhoE[%f,%f] = [%f,%f,%f,%f]\n',double(X),double(Y),double(subs(U,args,vals))); 
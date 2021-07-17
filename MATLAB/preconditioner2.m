syms gamma rho rhoU rhoV E;
U = [ rho, rhoU, rhoV, E];
q = (1/(2*rho))*(rhoU^2+rhoV^2);
p = (gamma-1)*(E-q);
s = log(p) - gamma*log(rho);
%s = log(p/(rho^gamma));
%S = -rho*s/(gamma-1);
S = -rho*s/(gamma-1);
W = [p,rhoU/rho,rhoV/rho,S];
dWdU = jacobian(W,U);
disp("dWdU Jacobian");
disp(dWdU);

%back substitute P
u=rhoU/rho;
v=rhoV/rho;
dWdUout=subs(dWdU,str2sym('E - (rhoU^2 + rhoV^2)/(2*rho)'),str2sym('P/(gamma-1)'));
syms u v gmg;
dWdUout=subs(dWdUout,rhoU/rho,u);
dWdUout=subs(dWdUout,rhoV/rho,v);
%dWdUout=subs(dWdUout,gamma/(gamma-1),gmg);

disp(dWdUout);
error("planned");

dUdW = inv(dWdU);
disp("dUdW Jacobian");
disp(dUdW);

syms alpha a BMr2Oa2 delta;
a = sqrt(gamma*p/rho);
P0 = [
    BMr2Oa2, 0, 0, -delta*BMr2Oa2;
    -alpha*rhoU/(rho^2*a^2), 1, 0, delta*alpha*rhoU/(rho^2*a^2);
    -alpha*rhoV/(rho^2*a^2), 0, 1, delta*alpha*rhoV/(rho^2*a^2);
    0, 0, 0, 1;
    ];
disp("pressure coordinates preconditioner");
P0 = simplify(P0);
disp(P0);

P0t = simplify(dUdW*P0*dWdU);

disp("transformed preconditioner in conservative variables");
disp(P0t);

disp("source code for transformed preconditioner");

ccode(P0t)
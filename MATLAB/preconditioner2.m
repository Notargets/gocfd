syms gamma rho rhoU rhoV E;
U = [ rho, rhoU, rhoV, E];
q = (1/(2*rho))*(rhoU^2+rhoV^2);
p = (gamma-1)*(E-q);
%S = log(p / (rho^gamma));
S = log(p/(rho^gamma));
%S = log(p) - gamma*log(rho);
W = [p,rhoU/rho,rhoV/rho,S];
dWdU = jacobian(W,U);
disp("dWdU Jacobian");
disp(dWdU);

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
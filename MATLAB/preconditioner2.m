syms gamma rho rhoU rhoV E a;
U = [ rho, rhoU, rhoV, E]; 
W = [ (gamma-1)*(E- (1/(2*rho))*(rhoU^2+rhoV^2)), rhoU/rho, rhoV/rho, ((gamma-1)/(rho^gamma))*(E- (1/(2*rho))*(rhoU^2+rhoV^2))];
dWdU = jacobian(W,U);
dUdW = inv(dWdU);
disp("dUdW Jacobian");
disp(dUdW);
syms alpha u v BMr2Oa2 delta;
P0 = [
    BMr2Oa2, 0, 0, -delta*BMr2Oa2;
    -alpha*u/(rho*a^2), 1, 0, delta*alpha*u/(rho*a^2);
    -alpha*v/(rho*a^2), 0, 1, delta*alpha*v/(rho*a^2);
    0, 0, 0, 1;
    ];
disp("pressure coordinates preconditioner");
disp(P0);

P0t = simplify(dUdW*P0*dWdU);

disp("transformed preconditioner in conservative variables");
disp(P0t);

disp("source code for transformed preconditioner");
ccode(P0t)
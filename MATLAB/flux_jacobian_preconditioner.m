% Check the jacobian function to ensure it's delivering proper answers for
% substitution variables
U = [ rho, rhoU, rhoV, E];
u = rhoU/rho;
v = rhoV/rho;
F = [ rhoU, rho*u^2+p, rho*u*v, u*(E+p) ];
G = [ rhoV, rho*u*v, rho*v^2+p, v*(E+p) ];
p=(gamma-1.)*(E-(rhoU^2+rhoV^2)/(2.*rho));

W = [ w0, w1, w2, w3 ];
pw=(gamma-1)*(w3-0.5*(w1^2+w2^2)/w0);
FW = [ w1, w1^2/w0+pw, w1*w2/w0, (w1/w0)*(w3+pw)];
dFWdW = simplify(jacobian(FW,W));
dFdU2 = subs(dFWdW,W,U);
dFdU = simplify(jacobian(F,U));

if simplify(dFdU2-dFdU) == 0 
    disp('correct answer')
end
dGdU = simplify(jacobian(G,U));

% Output the code to manifest each jacobian
ccode(dFdU)
ccode(dGdU)
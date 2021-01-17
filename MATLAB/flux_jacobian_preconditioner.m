% 2D Euler Equations - Preconditioner developement
%
% Theory: In short, we're taking a known RHS in our 2D Flux Reconstructed
% Galerkin Euler Solver and getting an algebraic expression for the
% derivative of that RHS relative to the change in the left hand side
% conservative variables. In this case:
%   dU/dt = -divergence([F(U),G(U)])
% So the RHS is just the divergence of the flux functions in the two
% directions X and Y. To do a Newton solver, we should get the Jacobian of
% the vector RHS(U) with respect to the vector U, where:
%   u = rhoU/rho, v = rhoV/rho % for clarity
%   U = [ rho, rhoU, rhoV, E ]
%   F = [ rhoU, rho*u^2+p, rhoU*v, u*(E+p) ]
%   G = [ rhoV, rhoV*u, rho*v^2+p, v*(E+p) ]
% so - the objective is to find:
%   dRHS/dU, which is a 4x4 matrix
% We assume that we can obtain the derivatives of the U field in X and Y,
% so, to compute the value of the Jacobian requires we compute Ux and Uy.
%
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
%ccode(dFdU)
%ccode(dGdU)

syms ux1 ux2 ux3 ux4 uy1 uy2 uy3 uy4;
dUdX = [ux1,ux2,ux3,ux4].';
dUdY = [uy1,uy2,uy3,uy4].';

DivFG = -(dFdU*dUdX+dGdU*dUdY);

P = jacobian(DivFG,U);
% Output the code to manifest the jacobian: DResidual/DU = P(U,Ux,Uy)
disp(ccode(simplify(P)));


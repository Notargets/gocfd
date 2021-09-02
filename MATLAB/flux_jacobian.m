% 2D Euler Equations - Preconditioner development
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
syms rho u v rhoU rhoV E gamma p;
p=(gamma-1.)*(E-(rhoU^2+rhoV^2)/(2.*rho));
U = [ rho, rhoU, rhoV, E];
F = [ rhoU, rhoU^2/rho+p, rhoU*rhoV/rho, (rhoU/rho)*(E+p) ];
G = [ rhoV, rhoV*rhoU/rho, rhoV^2/rho+p, (rhoV/rho)*(E+p) ];

dFdU = simplify(jacobian(F,U));
dGdU = simplify(jacobian(G,U));
syms u v;
dFdU = subs(dFdU,rhoU,str2sym('rho*u'));
dFdU = subs(dFdU,rhoV,str2sym('rho*v'));
dGdU = subs(dGdU,rhoU,str2sym('rho*u'));
dGdU = subs(dGdU,rhoV,str2sym('rho*v'));

% Output the code to manifest each jacobian
fprintf("dFdU - F component of flux jacobian\n%s\n",ccode(simplify(dFdU,"Steps",100)));
fprintf("dGdU - G component of flux jacobian\n%s\n",ccode(simplify(dGdU,"Steps",100)));

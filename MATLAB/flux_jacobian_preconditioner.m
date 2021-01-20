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
syms rho rhoU rhoV E gamma;
U = [ rho, rhoU, rhoV, E];
u = rhoU/rho;
v = rhoV/rho;
p=(gamma-1.)*(E-(rhoU^2+rhoV^2)/(2.*rho));
F = [ rhoU, rho*u^2+p, rho*u*v, u*(E+p) ];
G = [ rhoV, rho*u*v, rho*v^2+p, v*(E+p) ];

syms w0 w1 w2 w3;
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

% We need to calculate dRHS/dU, making the approximation that we can use
% the singular nodal values of U to compose dF/dU and dG/dU rather than sum
% all of the nodal calculations. The original system is this:
%       dU/dt(X_i) = - sum_over_I(Flux(X_i) * divergence(Psi(X_i))
% Flux in the above is a scalar DOF in an RT element, which for the
% solution points is a projection of the vector flux [F,G] onto X or Y
% directions, depening on the node, i.
%
% We are calculating a Jacobian of the RHS, dRHS/dU. Psi is a polynomial in
% X, and so doesn't depend on U, leaving the Flux, which is dependent on U.
% We want a value for dRHS/dU at location X_i, which would require we do
% this:
%       dRHS/dU = - sum_over_I(dFlux/dU(X_i) * divergence(Psi(X_i))
% But that would require Imax matrix summations per point within I, so
% instead we'll calculate this:
%       dRHS/dU ~= - dFlux/dU(X_i) * sum_over_I(divergence(Psi(X_i))
% The right most term is a scalar for each node i, so we'll be multiplying
% the local jacobian matrix by the local divergence, then calculating the
% inverse to compose:
%       dU ~= dt * [dRHS/dU]^-1 * RHS
%
% Note that the Flux_i term is a scalar projection of [F,G], so we need to
% compose d[F,G]/dU then project it onto the local flux directional vector
% for each internal point. For that, we need to calculate the directional
% vector, or obtain it from the calculations leading to the inital
% projection. We can calculate the flux direction using the projected flux
% from the two equations:
%       Flux_i = nx*F+ny*G, sqrt(nx^2+ny^2) = 1
%       ny = sqrt(1-nx^2) -> allows for ny +/-
% Note that we can't easily distinguish which direction of ny to use!
% Instead, we can obtain net individual flux directional values from the RT
% element interpolations, then we have:
%       dFlux/dU(X_i) = [dF/dU,dG/dU] dot [nx,ny]
% The dot product of [dF/dU,dG/dU] with [nx,ny} produces a single matrix
% containing the directional flux jacobian.
% We then need the inverse of the directional flux jacobian which we can
% then use as a preconditioner for the RHS as shown above.

syms nx ny;
Flux = dFdU*nx+dGdU*ny;
disp 'Inverse of directional flux: ';
[invFlux,sigma] = subexpr(simplify(inv(Flux)));
syms sigma2;
[invFlux2,sigma2] = subexpr(simplify(invFlux),sigma2);
%Output simplified preconditioner with subexpressions sigma and sigma2
disp 'preconditioner, ~= [dRHS/dU]^-1'
ccode(invFlux2)
disp 'sigma subexpression'
ccode(sigma)
disp 'sigma2 subexpression'
ccode(sigma2)
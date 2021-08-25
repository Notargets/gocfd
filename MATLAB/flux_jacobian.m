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
syms w1 w2;
dFluxDU = w1*dFdU+w2*dGdU;
dFluxDU = subs(dFluxDU,rhoU,u*rho);
dFluxDU = subs(dFluxDU,rhoV,v*rho);
dFluxDU = simplify(dFluxDU);
dFluxDU_inv = inv(dFluxDU);
fprintf("dFluxDU inverse\n[ %s , %s , %s , %s ]\n",dFluxDU_inv);

% Simplify using variable substitutions
syms x1 x2 x3 x4 x5 x6;
x1s = str2sym('(rho*gamma^2*u^3*w1^3 + rho*gamma^2*u^3*w1*w2^2 + rho*gamma^2*u^2*v*w1^2*w2 + rho*gamma^2*u^2*v*w2^3 + rho*gamma^2*u*v^2*w1^3 + rho*gamma^2*u*v^2*w1*w2^2 - 2*E*gamma^2*u*w1^3 - 2*E*gamma^2*u*w1*w2^2 + rho*gamma^2*v^3*w1^2*w2 + rho*gamma^2*v^3*w2^3 - 2*E*gamma^2*v*w1^2*w2 - 2*E*gamma^2*v*w2^3 - rho*gamma*u^3*w1^3 - rho*gamma*u^3*w1*w2^2 - rho*gamma*u^2*v*w1^2*w2 - rho*gamma*u^2*v*w2^3 - rho*gamma*u*v^2*w1^3 - rho*gamma*u*v^2*w1*w2^2 + 2*E*gamma*u*w1^3 + 2*E*gamma*u*w1*w2^2 - rho*gamma*v^3*w1^2*w2 - rho*gamma*v^3*w2^3 + 2*E*gamma*v*w1^2*w2 + 2*E*gamma*v*w2^3 + 2*rho*u^3*w1^3 + 6*rho*u^2*v*w1^2*w2 + 6*rho*u*v^2*w1*w2^2 + 2*rho*v^3*w2^3)');
x2s = str2sym('(u*w2 - v*w1 - gamma*u*w2 + gamma*v*w1)');
x3s = str2sym('rho*gamma^2*u^2 + rho*gamma^2*v^2 - 2*E*gamma^2 - 3*rho*gamma*u^2 + 2*rho*gamma*u*v - rho*gamma*v^2 + 2*E*gamma + 2*rho*u^2 + 2*rho*v^2');
x4s = str2sym('(rho*gamma^2*u^2*w1^2 + rho*gamma^2*u^2*w2^2 + rho*gamma^2*v^2*w1^2 + rho*gamma^2*v^2*w2^2 - 2*E*gamma^2*w1^2 - 2*E*gamma^2*w2^2 + 2*E*gamma*w1^2 + 2*E*gamma*w2^2 + 3*rho*u^2*w1^2 - rho*u^2*w2^2 + 8*rho*u*v*w1*w2 - rho*v^2*w1^2 + 3*rho*v^2*w2^2)');
x5s = str2sym('rho*gamma^2*u^3 + rho*gamma^2*u^2*v + rho*gamma^2*u*v^2 - 2*E*gamma^2*u + rho*gamma^2*v^3 - 2*E*gamma^2*v + 2*E*gamma*u + 2*E*gamma*v + rho*u^3 + rho*v^3');
d2 = subs(dFluxDU_inv, x1s, x1);
d3 = subs(d2,x2s,x2);
d4 = subs(d3,x3s,x3);
d5 = subs(d4,x4s,x4);
d6 = subs(d5,x5s,x5);
%d2 = subs(dFluxDU_inv,str2sym('(rho*gamma^2*u^3 + rho*gamma^2*u^2*v + rho*gamma^2*u*v^2 - 2*E*gamma^2*u + rho*gamma^2*v^3 - 2*E*gamma^2*v - rho*gamma*u^3 - rho*gamma*u^2*v - rho*gamma*u*v^2 + 2*E*gamma*u - rho*gamma*v^3 + 2*E*gamma*v + rho*u^3 + 3*rho*u^2*v + 3*rho*u*v^2 + rho*v^3)'),x1);
%d3 = subs(d2,str2sym('rho*gamma^2*u^2 + rho*gamma^2*v^2 - 2*E*gamma^2 - rho*gamma*u^2 + 2*rho*gamma*u*v - 3*rho*gamma*v^2 + 2*E*gamma + 2*rho*u^2 + 2*rho*v^2'),x2);
%d4 = subs(d3,str2sym('rho*gamma^2*u^2 + rho*gamma^2*v^2 - 2*E*gamma^2 - 3*rho*gamma*u^2 + 2*rho*gamma*u*v - rho*gamma*v^2 + 2*E*gamma + 2*rho*u^2 + 2*rho*v^2'),x3);
%d5 = subs(d4,str2sym('rho*gamma^2*u^2 + rho*gamma^2*v^2 - 2*E*gamma^2 + 2*E*gamma + rho*u^2 + 4*rho*u*v + rho*v^2'),x4);
%d6 = subs(d5,str2sym('rho*gamma^2*u^3 + rho*gamma^2*u^2*v + rho*gamma^2*u*v^2 - 2*E*gamma^2*u + rho*gamma^2*v^3 - 2*E*gamma^2*v + 2*E*gamma*u + 2*E*gamma*v + rho*u^3 + rho*v^3'),x5);
dFluxDU_inv = simplify(d6);

% Output the code to manifest each jacobian
fprintf("dFdU - F component of flux jacobian\n%s\n",ccode(simplify(dFdU)));
fprintf("dGdU - G component of flux jacobian\n%s\n",ccode(simplify(dGdU)));
fprintf("dFluxdU - dFdU+dGdU sum of flux jacobians projected onto [w1,w2]\n%s\n",ccode(simplify(dFluxDU)));
fprintf("dFluxdU_Inv - Inverse\n%s\n",ccode(x1s),ccode(x2s),ccode(x3s),ccode(x4s),ccode(x5s),ccode(simplify(dFluxDU_inv)));

x1s = simplify(x1s);
x2s = simplify(x2s);
x3s = simplify(x3s);
x4s = simplify(x4s);
x5s = simplify(x5s);

fprintf("dFluxDU Inverse\n[ %s , %s , %s , %s ]\n",dFluxDU_inv);
fprintf("x1\n%s\n",x1s);
fprintf("x2\n%s\n",x2s);
fprintf("x3\n%s\n",x3s);
fprintf("x4\n%s\n",x4s);
fprintf("x5\n%s\n",x5s);

%syms u v;
%D = eig(simplify(dFluxDU));
%D = subs(D,rhoU,u*rho);
%D = subs(D,rhoV,v*rho);
%fprintf("Eigenvalues of dFluxdU - w1 and w2 are the nodal basis unit vector values\n%s\n",ccode(simplify(D)));

% We need to calculate dRHS/dU, making the approximation that we can use
% the singular nodal values of U to compose dF/dU and dG/dU rather than sum
% all of the nodal calculations. The original system is this:
%       dU/dt(X_i) = - sum_over_I(Flux(X_i) * divergence(Psi(X_i))
% Flux in the above is a scalar DOF in an RT element, which for the
% interior solution points is the magnitude of the vector flux [F,G].
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
% We then need the inverse of the directional flux jacobian which we can
% then use as a preconditioner for the RHS as shown above.
%Flux = sqrt(F.^2+G.^2);
%dFluxdU = jacobian(simplify(Flux),U);
% To ease the use of the preconditioner, we'll divide by the Flux so that
% we can compose the dRHS/dU as:
%    dRHS/dU ~= (dFluxdU/Flux)*(sum(Flux*divergence(Psi))
% Which is approximately equal to:
%    dRHS/dU ~= dFluxdU*sum(divergence(Psi))
%precon = simplify(dFluxDU);
% Reduce and simplify the result, output code
%for m = 1:19
%    name = ['Ssig' num2str(m,'%d')];
%    name2 = name;
%    sym name;
%    [precon,name] = subexpr(simplify(precon),name);
%    fprintf ("%s = %s\n", char(name2), ccode(name));
%end
%fprintf("%s\n",ccode(precon));


% Instead of below, we'll do a numerical inverse - it will be cheaper
% disp 'Inverse of directional flux: ';
% [invFlux,sigma] = subexpr(simplify(inv(dFluxdU)));
% syms sigma2;
% [invFlux2,sigma2] = subexpr(simplify(invFlux),sigma2);
% %Output simplified preconditioner with subexpressions sigma and sigma2
% disp 'preconditioner, ~= [dRHS/dU]^-1'
% ccode(invFlux2)
% disp 'sigma subexpression'
% ccode(sigma)
% disp 'sigma2 subexpression'
% ccode(sigma2)
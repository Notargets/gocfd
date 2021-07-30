syms gamma rho rhoU rhoV E;
U = [ rho, rhoU, rhoV, E];
q = (1/(2*rho))*(rhoU^2+rhoV^2);
p = (gamma-1)*(E-q);
s = log(p) - gamma*log(rho);
%s = log(p/(rho^gamma));
%S = -rho*s/(gamma-1);
S = s;
W = [p,rhoU/rho,rhoV/rho,S];
dWdU = jacobian(W,U);
%disp("dWdU Jacobian");
%disp(dWdU);

%back substitute P
syms u v gmg P;
u=rhoU/rho;
v=rhoV/rho;
dWdUout=subs(dWdU,str2sym('E - (rhoU^2 + rhoV^2)/(2*rho)'),str2sym('P/(gamma-1)'));
dWdUout=subs(dWdUout,rhoU/rho,u);
dWdUout=subs(dWdUout,rhoV/rho,v);
%dWdUout=subs(dWdUout,gamma/(gamma-1),gmg);
dWdU = dWdUout;
%disp(dWdU);
dUdW = inv(dWdUout);
%disp("dUdW Jacobian");
%disp(dUdW);

syms alpha a B2 delta;
a = sqrt(gamma*p/rho);
%P0 = [
%    BMr2Oa2, 0, 0, -delta*BMr2Oa2;
%    -alpha*rhoU/(rho^2*a^2), 1, 0, delta*alpha*rhoU/(rho^2*a^2);
%    -alpha*rhoV/(rho^2*a^2), 0, 1, delta*alpha*rhoV/(rho^2*a^2);
%    0, 0, 0, 1;
%    ];
syms u v;
rhoU = u*rho;
rhoV = v*rho;
P0 = [
    B2, 0, 0, 0;
    0, 1, 0, 0;
    0, 0, 1, 0;
    0, 0, 0, 1;
    ];
disp("pressure coordinates preconditioner");
%P0 = simplify(P0);
disp(P0);

dWdU = subs(dWdU,str2sym('rhoU'),rho*u);
dWdU = subs(dWdU,str2sym('rhoV'),rho*v);
dUdW = subs(dUdW,str2sym('rhoU'),rho*u);
dUdW = subs(dUdW,str2sym('rhoV'),rho*v);
%q = (1/2)*rho*qq;
dUdW = simplify(dUdW,'Steps',100);
dWdU = simplify(dWdU,'Steps',100);
dUdW = subs(dUdW,str2sym('u^2+v^2'),str2sym('qq'));
dWdU = subs(dWdU,str2sym('u^2+v^2'),str2sym('qq'));
%qqq = qq*(gamma-1)/2
syms c2 qqq;
dUdW = subs(dUdW,str2sym('(qq*(gamma-1))/2'),qqq);
dWdU = subs(dWdU,str2sym('(qq*(gamma-1))/2'),qqq);
dUdW = subs(dUdW,str2sym('rho/(P*gamma)'),str2sym('1/c2'));
dUdW = subs(dUdW,str2sym('rho/(gamma)'),str2sym('1/(c2*P)'));
dUdW = subs(dUdW,str2sym('qq/c2'),str2sym('m2'));
%dUdW(4,1) verified by hand to be equal to h/c2 where h=c2/(gamma-1)+qq/2
syms h;
dUdW(4,1) = h/c2;
%dWdU(4,1) verified by hand to be equal to (qqq-c2)/P
dWdU(4,1) = (qqq-c2)/P;
disp('dUdW');
disp(dUdW);
disp('dWdU');
disp(dWdU);
disp('dWdU with pressure error ommitted (like the Turkel paper)');
dWdU(4,:) =  dWdU(4,:)*P;
disp(dWdU);
%error('planned');
P0t = simplify(dUdW*P0*dWdU,'Steps',100);
P0t = subs(P0t,str2sym('rho/(P*gamma)'),str2sym('1/c2'));
P0t = subs(P0t,str2sym('rho/(gamma)'),str2sym('1/(c2*P)'));
P0t = subs(P0t,str2sym('qq/c2'),str2sym('m2'));
syms q;
P0t = subs(P0t,str2sym('(rho*u^2)/2 + (rho*v^2)/2'),q);
%q=(rho*u^2)/2 + (rho*v^2)/2;
P0t = simplify(P0t,'Steps',100);
%syms qq;
%P0t = subs(P0t,str2sym('(u^2 + v^2)'),qq);
%qq = (u^2 + v^2);
%P0t = simplify(P0t,'Steps',100);
disp("transformed preconditioner in conservative variables");
disp(P0t);

disp("source code for transformed preconditioner");

ccode(P0t)
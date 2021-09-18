flow_variables

syms rhoL rhoR h0L h0R uL uR vL vR wL wR;

uRL = roeAvg(uL,uR);
vRL = roeAvg(vL,vR);
wRL = roeAvg(wL,wR);
h0RL = roeAvg(h0L,h0R);

rhoRL = sqrt(rhoL*rhoR);

magURL = str2sym('sqrt(uRL^2+vRL^2+wRL^2)');
cRL = str2sym('sqrt((Gamma-1)*(h0RL-magURL^2/2))');


disp("loaded roe averaged variables and function roeAvg(sL,sR)");

function p = roeAvg(sL,sR)
% This is used to calculate the three velocity components and the enthalpy
    syms rhoL rhoR;
    p = (sL*sqrt(rhoL) + sR*sqrt(rhoR))/(sqrt(rhoL)+sqrt(rhoR));
end
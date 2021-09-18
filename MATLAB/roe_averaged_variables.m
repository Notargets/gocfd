flow_variables

syms rhoL rhoR h0L h0R uL uR;

uRL = roeAvg(uL,uR);
h0RL = roeAvg(h0L,h0R);

rhoRL = sqrt(rhoL*rhoR);

cRL = sqrt((Gamma-1)*(h0RL-uRL^2/2));

disp("loaded roe averaged variables and function roeAvg(sL,sR)");

function p = roeAvg(sL,sR)
    syms rhoL rhoR;
    p = (sL*sqrt(rhoL) + sR*sqrt(rhoR))/(sqrt(rhoL)+sqrt(rhoR));
end
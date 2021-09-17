syms rho u e0 E0 T Cp Cv Gamma M;
U = [rho rho*u rho*e0];
U = [rho rho*u E0];
% Perfect Gas
R = Cp - Cv;
Gamma = Cp/Cv;
P = rho*R*T;
% Specific Energy
e = e0 - u^2/2;
e = Cv*T;
% Pressure
q = (1/2)*rho*u^2;
P = (Gamma-1)*(rho*e0-q);
% Stagnation Enthalpy
h0 = e0 + P/rho;
% Specific Enthalpy
h = h0 - q/rho;
% Stagnation Pressure
P0 = P*((1+(Gamma-1)/2)*M^2)^(Gamma/(Gamma-1));
% Stagnation Temperature
T0 = T + (q/(rho*Cp));
% Speed of Sound
c = sqrt(Gamma*R*T);
% Mach Number
M = u/c;
% Enthalpy and Total Enthalpy
H = rho*h;
H0 = rho*h0;
% Energy and Total Energy
E = rho*e;
E0 = rho*e0;
% Note that the symbol "E" is commonly used to refer to Total Energy when
% seen in conservative variables
disp("loaded flow variables");

syms rho u e0 E0 T Cp Cv Gamma M;
U = [rho rho*u rho*e0];
U = [rho rho*u E0];
% Perfect Gas
% R = Cp - Cv;
R = 1;
%Gamma = Cp/Cv;
P = rho*R*T;
% Specific Energy
e = Cv*T;
e = str2sym('e0 - u^2/2');
% Pressure
P = str2sym('(E0 - q)*(Gamma - 1)');
q = (rho*u^2)/2;
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
c = str2sym('sqrt(Gamma*P/rho)');
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

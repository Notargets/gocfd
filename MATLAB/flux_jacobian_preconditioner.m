p=(gamma-1.)*(E-(rhoU^2+rhoV^2)/(2.*rho));
U = [ rho, rhoU, rhoV, E];
F = [ rhoU, rho*u^2+p, rho*u*v, u*(E+p) ];
G = [ rhoV, rho*u*v, rho*v^2+p, v*(E+p) ];

W = [ w0, w1, w2, w3 ];
pw=(gamma-1)*(w3-0.5*(w1^2+w2^2)/w0);

FW = [ w1, w1^2/w0+pw, w2^2/w0+pw, (w2/w0)*(w3+p)];
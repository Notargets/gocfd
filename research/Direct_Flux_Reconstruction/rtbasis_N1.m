%-- 8/11/2020 10:40 AM --%
syms a1 a2 a3 a4 a5 a6 a7 a8
syms x y
psi1 = [ a1+a2*x+a3*y+a7*x^2+a8*x*y, a4+a5*x+a6*y+a7*x*y+a8*y^2];
n1=[1,0];
n2=[0,1];
n3 = [ 1/sqrt(2), 1/sqrt(2) ];
n4 = n3;
n5=[-1,0];
n6 = n5;
n7 = [0,-1];
n8=n7;

aa1=-0.33333333;
aa2=0.44721360;

%Interior
r1=[aa1,aa1];
r2=r1;
%Edge 1
r3=[aa2,-aa2];
r4=[-aa2,aa2];
%Edge 2
r5=[-1,-aa2];
r6=[-1,aa2];
%Edge 3
r7=[-aa2,-1];
r8=[aa2,-1];

eq1=subs(psi1,{x,y}, r1)*n1';
eq2=subs(psi1,{x,y}, r2)*n2';
eq3=subs(psi1,{x,y}, r3)*n3';
eq4=subs(psi1,{x,y}, r4)*n4';
eq5=subs(psi1,{x,y}, r5)*n5';
eq6=subs(psi1,{x,y}, r6)*n6';
eq7=subs(psi1,{x,y}, r7)*n7';
eq8=subs(psi1,{x,y}, r8)*n8';
A = equationsToMatrix([eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8]);
digits(4);
disp("A =");
disp(vpa(A));
Ainv = inv(A);
disp("Ainv = ");
disp(vpa(inv(A)));

%We have coefficients for the basis, extract the basis functions 1-8
b1 = subs(psi1, {a1,a2,a3,a4,a5,a6,a7,a8}, Ainv(:,1)');
b2 = subs(psi1, {a1,a2,a3,a4,a5,a6,a7,a8}, Ainv(:,2)');
b3 = subs(psi1, {a1,a2,a3,a4,a5,a6,a7,a8}, Ainv(:,3)');
b4 = subs(psi1, {a1,a2,a3,a4,a5,a6,a7,a8}, Ainv(:,4)');
b5 = subs(psi1, {a1,a2,a3,a4,a5,a6,a7,a8}, Ainv(:,5)');
b6 = subs(psi1, {a1,a2,a3,a4,a5,a6,a7,a8}, Ainv(:,6)');
b7 = subs(psi1, {a1,a2,a3,a4,a5,a6,a7,a8}, Ainv(:,7)');
b8 = subs(psi1, {a1,a2,a3,a4,a5,a6,a7,a8}, Ainv(:,8)');

%Construct the Vandermonde Matrices for each of r and s directions
xx = [[r1]; [r2]; [r3]; [r4]; [r5]; [r6]; [r7]; [r8];];
b1eval = evalb(b1, x, y, xx);
b2eval = evalb(b2, x, y, xx);
b3eval = evalb(b3, x, y, xx);
b4eval = evalb(b4, x, y, xx);
b5eval = evalb(b5, x, y, xx);
b6eval = evalb(b6, x, y, xx);
b7eval = evalb(b7, x, y, xx);
b8eval = evalb(b8, x, y, xx);

V1 = [b1eval(:,1), b2eval(:,1), b3eval(:,1), b4eval(:,1), b5eval(:,1), b6eval(:,1), b7eval(:,1), b8eval(:,1)];
V2 = [b1eval(:,2), b2eval(:,2), b3eval(:,2), b4eval(:,2), b5eval(:,2), b6eval(:,2), b7eval(:,2), b8eval(:,2)];

%Compute derivative matrices, Dr and Ds
%Re-use all the above by redefining the basis function as the derivative
%(Dr first):
psi1 = [ a1+a2*x+a3*y+a7*x^2+a8*x*y, a4+a5*x+a6*y+a7*x*y+a8*y^2];

psi1r = [ a2+2*a7*x+a8*y, a5+a7*y];
b1 = subs(psi1r, {a1,a2,a3,a4,a5,a6,a7,a8}, Ainv(:,1)');
b2 = subs(psi1r, {a1,a2,a3,a4,a5,a6,a7,a8}, Ainv(:,2)');
b3 = subs(psi1r, {a1,a2,a3,a4,a5,a6,a7,a8}, Ainv(:,3)');
b4 = subs(psi1r, {a1,a2,a3,a4,a5,a6,a7,a8}, Ainv(:,4)');
b5 = subs(psi1r, {a1,a2,a3,a4,a5,a6,a7,a8}, Ainv(:,5)');
b6 = subs(psi1r, {a1,a2,a3,a4,a5,a6,a7,a8}, Ainv(:,6)');
b7 = subs(psi1r, {a1,a2,a3,a4,a5,a6,a7,a8}, Ainv(:,7)');
b8 = subs(psi1r, {a1,a2,a3,a4,a5,a6,a7,a8}, Ainv(:,8)');

%Construct the Vandermonde Matrices for each of r and s directions
b1eval = evalb(b1, x, y, xx);
b2eval = evalb(b2, x, y, xx);
b3eval = evalb(b3, x, y, xx);
b4eval = evalb(b4, x, y, xx);
b5eval = evalb(b5, x, y, xx);
b6eval = evalb(b6, x, y, xx);
b7eval = evalb(b7, x, y, xx);
b8eval = evalb(b8, x, y, xx);

Dr1 = [b1eval(:,1), b2eval(:,1), b3eval(:,1), b4eval(:,1), b5eval(:,1), b6eval(:,1), b7eval(:,1), b8eval(:,1)];
Dr2 = [b1eval(:,2), b2eval(:,2), b3eval(:,2), b4eval(:,2), b5eval(:,2), b6eval(:,2), b7eval(:,2), b8eval(:,2)];

psi1 = [ a1+a2*x+a3*y+a7*x^2+a8*x*y, a4+a5*x+a6*y+a7*x*y+a8*y^2];

psi1s = [ a3+a8*x, a6+a7*x+2*a8*y];
b1 = subs(psi1s, {a1,a2,a3,a4,a5,a6,a7,a8}, Ainv(:,1)');
b2 = subs(psi1s, {a1,a2,a3,a4,a5,a6,a7,a8}, Ainv(:,2)');
b3 = subs(psi1s, {a1,a2,a3,a4,a5,a6,a7,a8}, Ainv(:,3)');
b4 = subs(psi1s, {a1,a2,a3,a4,a5,a6,a7,a8}, Ainv(:,4)');
b5 = subs(psi1s, {a1,a2,a3,a4,a5,a6,a7,a8}, Ainv(:,5)');
b6 = subs(psi1s, {a1,a2,a3,a4,a5,a6,a7,a8}, Ainv(:,6)');
b7 = subs(psi1s, {a1,a2,a3,a4,a5,a6,a7,a8}, Ainv(:,7)');
b8 = subs(psi1s, {a1,a2,a3,a4,a5,a6,a7,a8}, Ainv(:,8)');

%Construct the Vandermonde Matrices for each of r and s directions
b1eval = evalb(b1, x, y, xx);
b2eval = evalb(b2, x, y, xx);
b3eval = evalb(b3, x, y, xx);
b4eval = evalb(b4, x, y, xx);
b5eval = evalb(b5, x, y, xx);
b6eval = evalb(b6, x, y, xx);
b7eval = evalb(b7, x, y, xx);
b8eval = evalb(b8, x, y, xx);

Ds1 = [b1eval(:,1), b2eval(:,1), b3eval(:,1), b4eval(:,1), b5eval(:,1), b6eval(:,1), b7eval(:,1), b8eval(:,1)];
Ds2 = [b1eval(:,2), b2eval(:,2), b3eval(:,2), b4eval(:,2), b5eval(:,2), b6eval(:,2), b7eval(:,2), b8eval(:,2)];

function bval = evalb(basis, x, y, xx)
    b = basis;
    bval = [
        subs(b, {x,y}, xx(1,:))
        subs(b, {x,y}, xx(2,:))
        subs(b, {x,y}, xx(3,:))
        subs(b, {x,y}, xx(4,:))
        subs(b, {x,y}, xx(5,:))
        subs(b, {x,y}, xx(6,:))
        subs(b, {x,y}, xx(7,:))
        subs(b, {x,y}, xx(8,:))
    ];
end


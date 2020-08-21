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
b = b1;
b1eval = [
    subs(b, {x,y}, r1),
    subs(b, {x,y}, r2),
    subs(b, {x,y}, r3),
    subs(b, {x,y}, r4),
    subs(b, {x,y}, r5),
    subs(b, {x,y}, r6),
    subs(b, {x,y}, r7),
    subs(b, {x,y}, r8)
];
b = b2;
b2eval = [
    subs(b, {x,y}, r1),
    subs(b, {x,y}, r2),
    subs(b, {x,y}, r3),
    subs(b, {x,y}, r4),
    subs(b, {x,y}, r5),
    subs(b, {x,y}, r6),
    subs(b, {x,y}, r7),
    subs(b, {x,y}, r8)
];
b = b3;
b3eval = [
    subs(b, {x,y}, r1),
    subs(b, {x,y}, r2),
    subs(b, {x,y}, r3),
    subs(b, {x,y}, r4),
    subs(b, {x,y}, r5),
    subs(b, {x,y}, r6),
    subs(b, {x,y}, r7),
    subs(b, {x,y}, r8)
];
b = b4;
b4eval = [
    subs(b, {x,y}, r1),
    subs(b, {x,y}, r2),
    subs(b, {x,y}, r3),
    subs(b, {x,y}, r4),
    subs(b, {x,y}, r5),
    subs(b, {x,y}, r6),
    subs(b, {x,y}, r7),
    subs(b, {x,y}, r8)
];
b = b5;
b5eval = [
    subs(b, {x,y}, r1),
    subs(b, {x,y}, r2),
    subs(b, {x,y}, r3),
    subs(b, {x,y}, r4),
    subs(b, {x,y}, r5),
    subs(b, {x,y}, r6),
    subs(b, {x,y}, r7),
    subs(b, {x,y}, r8)

];
b = b6;
b6eval = [
    subs(b, {x,y}, r1),
    subs(b, {x,y}, r2),
    subs(b, {x,y}, r3),
    subs(b, {x,y}, r4),
    subs(b, {x,y}, r5),
    subs(b, {x,y}, r6),
    subs(b, {x,y}, r7),
    subs(b, {x,y}, r8)

];
b = b7;
b7eval = [
    subs(b, {x,y}, r1),
    subs(b, {x,y}, r2),
    subs(b, {x,y}, r3),
    subs(b, {x,y}, r4),
    subs(b, {x,y}, r5),
    subs(b, {x,y}, r6),
    subs(b, {x,y}, r7),
    subs(b, {x,y}, r8)
];
b = b8;
b8eval = [
    subs(b, {x,y}, r1),
    subs(b, {x,y}, r2),
    subs(b, {x,y}, r3),
    subs(b, {x,y}, r4),
    subs(b, {x,y}, r5),
    subs(b, {x,y}, r6),
    subs(b, {x,y}, r7),
    subs(b, {x,y}, r8)
];

V1 = [b1eval(:,1), b2eval(:,1), b3eval(:,1), b4eval(:,1), b5eval(:,1), b6eval(:,1), b7eval(:,1), b8eval(:,1)];
V2 = [b1eval(:,2), b2eval(:,2), b3eval(:,2), b4eval(:,2), b5eval(:,2), b6eval(:,2), b7eval(:,2), b8eval(:,2)];

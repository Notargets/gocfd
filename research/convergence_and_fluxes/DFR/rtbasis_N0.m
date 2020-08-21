%-- 8/11/2020 10:40 AM --%
syms a1 a2 a3
syms x y
psi1 = [ a1+a3*x, a2+a3*y]
n1 = [ 1/sqrt(2), 1/sqrt(2) ]
n2 =[-1,0]
n3 = [0,-1]

%Edge 1
r1=[0,0]
%Edge 2
r2=[-1,0]
%Edge 3
r3=[0,-1]

eq1=subs(psi1,{x,y}, r1)*n1'
eq2=subs(psi1,{x,y}, r2)*n2'
eq3=subs(psi1,{x,y}, r3)*n3'
A = equationsToMatrix([eq1,eq2,eq3])
digits(4)
disp("A =")
disp(vpa(A))
disp("Ainv = ")
disp(vpa(inv(A)))
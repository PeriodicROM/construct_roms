function f = hk4_rhs(x,a,s,k,R) 
% Vars: x = [p0,1 p1,1 t0,2 t1,1] 
 
r = R^(1/2);

f = [-x(1)*s/r,...
 -x(2)*s*(k^2 + 1)/r + a(4)*k*s*x(4)/(a(2)*(k^2 + 1)),...
a(2)*a(4)*k*x(2)*x(4)/(2*a(3)) - 4*x(3)/r,...
 -a(2)*a(3)*k*x(2)*x(3)/a(4) + a(2)*k*x(2)/a(4) + x(4)*(-k^2 - 1)/r];
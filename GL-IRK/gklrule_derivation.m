clear
close all
clc

%%

syms t x0 w0 x1 w1 x2 w2 x3 w3 x4 w4 real

x = [x0 x1 x2 x3 x4];
w = [w0 w1 w2 w3 w4];

% n = 1 (1 point gauss, degree 2n-1 = 1)
% add three points (degree 3n+1 = 4)

xg = [1/2-sqrt(sym(3))/6, 1/2+sqrt(sym(3))/6];
wg = [sym(1)/2, sym(1)/2];

eqn = [dot(w, x.^0) == int(t.^0, 0, 1); ...
       dot(w, x.^1) == int(t.^1, 0, 1); ...
       dot(w, x.^2) == int(t.^2, 0, 1); ...
       dot(w, x.^3) == int(t.^3, 0, 1); ...
       dot(w, x.^4) == int(t.^4, 0, 1); ...
       dot(w, x.^5) == int(t.^5, 0, 1); ...
       dot(w, x.^6) == int(t.^6, 0, 1); ...
       dot(w, x.^7) == int(t.^7, 0, 1); ...
       x1 == xg(1); ...
       x3 == xg(2)];

assume(x0<x1<x2<x3<x4)
soln = solve(eqn, [x w]);

x0 = simplify(soln.x0);
x1 = simplify(soln.x1);
x2 = simplify(soln.x2);
x3 = simplify(soln.x3);
x4 = simplify(soln.x4);

x = [x0 x1 x2 x3 x4];

w0 = simplify(soln.w0);
w1 = simplify(soln.w1);
w2 = simplify(soln.w2);
w3 = simplify(soln.w3);
w4 = simplify(soln.w4);

w = [w0 w1 w2 w3 w4];

%%


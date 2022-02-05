clear
close all
clc

%%

syms t x0 w0 x1 w1 x2 w2 x3 w3 x4 w4 real

x = [x0 x1 x2];
w = [w0 w1 w2];

eqn = [dot(w, x.^0) == int(int(t.^0, 0, t), 0, 1); ...
       dot(w, x.^1) == int(int(t.^1, 0, t), 0, 1); ...
       dot(w, x.^2) == int(int(t.^2, 0, t), 0, 1); ...
       dot(w, x.^3) == int(int(t.^3, 0, t), 0, 1); ...
       dot(w, x.^4) == int(int(t.^4, 0, t), 0, 1); ...
       dot(w, x.^5) == int(int(t.^5, 0, t), 0, 1)];

% assume(x0<x1<x2)
soln = solve(eqn, [x w]);

x0 = simplify(soln.x0(1));
x1 = simplify(soln.x1(1));
x2 = simplify(soln.x2(1));
% x3 = simplify(soln.x3(1));
% x4 = simplify(soln.x4);

x = [x0 x1 x2];
[x, i] = sort(x);

w0 = simplify(soln.w0(1));
w1 = simplify(soln.w1(1));
w2 = simplify(soln.w2(1));
% w3 = simplify(soln.w3(1));
% w4 = simplify(soln.w4);

w = [w0 w1 w2];
w = w(i);

%%


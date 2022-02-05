clear
close all
clc

%%

syms t x0 w0 x1 w1 x2 w2 x3 w3 x4 w4 real

x = [x0 x1];
w = [w0 w1];

% eqn = [dot(w, x.^0) == int((t.^0)*sqrt((1-t)/(1+t))*2/pi, -1, 1); ...
%        dot(w, x.^1) == int((t.^1)*sqrt((1-t)/(1+t))*2/pi, -1, 1); ...
%        dot(w, x.^2) == int((t.^2)*sqrt((1-t)/(1+t))*2/pi, -1, 1); ...
%        dot(w, x.^3) == int((t.^3)*sqrt((1-t)/(1+t))*2/pi, -1, 1); ...
%        dot(w, x.^4) == int((t.^4)*sqrt((1-t)/(1+t))*2/pi, -1, 1); ...
%        dot(w, x.^5) == int((t.^5)*sqrt((1-t)/(1+t))*2/pi, -1, 1); ...
%        dot(w, x.^6) == int((t.^6)*sqrt((1-t)/(1+t))*2/pi, -1, 1); ...
%        dot(w, x.^7) == int((t.^7)*sqrt((1-t)/(1+t))*2/pi, -1, 1); ...
%        dot(w, x.^8) == int((t.^8)*sqrt((1-t)/(1+t))*2/pi, -1, 1); ...
%        dot(w, x.^9) == int((t.^9)*sqrt((1-t)/(1+t))*2/pi, -1, 1)];

eqn = [dot(w, x.^0) == int((t.^0)*sqrt((1-t)/(1+t))*2/pi, -1, 1); ...
       dot(w, x.^1) == int((t.^1)*sqrt((1-t)/(1+t))*2/pi, -1, 1); ...
       dot(w, x.^2) == int((t.^2)*sqrt((1-t)/(1+t))*2/pi, -1, 1); ...
       dot(w, x.^3) == int((t.^3)*sqrt((1-t)/(1+t))*2/pi, -1, 1)];

% assume(x0<x1<x2<x3<x4)
soln = solve(eqn, [x w]);

x0 = simplify(soln.x0(1));
x1 = simplify(soln.x1(1));
% x2 = simplify(soln.x2(1));
% x3 = simplify(soln.x3(1));
% x4 = simplify(soln.x4(1));

x = [x0 x1];

w0 = simplify(soln.w0(1));
w1 = simplify(soln.w1(1));
% w2 = simplify(soln.w2(1));
% w3 = simplify(soln.w3(1));
% w4 = simplify(soln.w4(1));

w = [w0 w1];

%%


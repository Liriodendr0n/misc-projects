clear
close all
clc

%%

syms x h real

syms f(x)

c = 2;

s = sym([-c, -1, 0, 1, c]);
% s = sym([-1, -c, 0, c, 1]);

J = vander(s);

a = flip(J'*J\J', 1);

% d0 = a(1,:);
d1 = a(2,:);
% d2 = a(3,:);
% d3 = a(4,:);
% 

d(h) = range(s)*d1*f(h*s')/(h);

order = 6;

T = taylor(d(h), h, 'order', order)
% Tf = taylor(f(h), h, 'order', order)
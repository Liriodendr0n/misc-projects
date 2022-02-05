clear
close all
clc

%% polynomial downwash

% g(t) = 1/pi * sqrt((1-x)/(1+x)) * x^N
% w(t) = p(t)

syms t

N = 2;

f(t) = 1/sym(pi) * sqrt((1-t)/(1+t));

nCk = @(n, k) factorial(n)./(factorial(k).*factorial(n-k));


m = sym(0:N);

pN = nCk(2*floor(m/2), floor(m/2))./(4.^floor(m/2) .* (-1).^(m+1));

p(t) = dot(pN, t.^(N:-1:0));

% figure
% hold on
% fplot(f, [-1, 1])
% fplot(p, [-1, 1])

%% quadrature rule

syms x0 xi0 w0 x1 xi1 w1



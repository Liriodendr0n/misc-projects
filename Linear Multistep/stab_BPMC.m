clear
close all
clc

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

syms x w theta
%% Bashforth Predictor - Moulton Corrector
%% find coefficients

% stencil = 0;
% stencil = [0 1];
% stencil = [0 1 2];
% stencil = [0 1 2 3];
% stencil = [0 1 3 5 6];
% stencil = [0 1 3 5 7 8];
% stencil = [0 1 3 5 7 9 10];
% stencil = [0 1 3 6 9 12 14 15];
% stencil = [0 1 3 6 9 12 15 17 18];
% stencil = [0 1 3 6 9 12 15 18 20 21];
% stencil = [0 1 3 6 9 13 17 20 23 25 26];
% stencil = [0 1 3 6 9 13 17 21 24 27 29 30];
stencil = [0 1 3 6 10 14 18 22 26 29 31 32];
% stencil = [0 1 3 6 10 14 18 22 26 30 33 35 36];
% stencil = [0 1 3 6 10 14 19 24 29 33 37 40 42 43];
% stencil = [0 1 3 6 10 14 19 24 29 34 38 42 45 47 48];
% stencil = [0 1 3 6 10 14 19 24 29 34 39 43 47 50 52 53];

n = length(stencil);

lagrpol = @(j) prod((x - stencil(1:end ~= j)')./(stencil(j) - stencil(1:end ~= j)'), 1);
for i = 1:n
    if n == 1
        l = 1;
    else
        l(i) = lagrpol(n-i+1);
    end
end

alpha = [zeros(1, max(stencil)), -1, 1];
betai = int(l, x, -1, 0);
gamma = [zeros(1, max(stencil)), -1, 1];
deltai = int(l, x, 0, 1);

for i = 1:length(stencil)
    beta(stencil(i)+1) = betai(i);
    delta(stencil(i)+1) = deltai(i);
end
beta = [beta 0];
delta = [0 delta];

k = length(beta)-1;

s(theta) = exp(1i*theta);

A(w) = delta(end).*sum(beta(1:end-1).*w.^(0:k-1));
B(w) = sum((delta(1:end-1)-delta(end)*alpha(1:end-1)).*w.^(0:k-1));
C(w) = -sum(gamma.*w.^(0:k));

h1 = (-B + sqrt(B.^2 - 4.*A.*C))./(2.*A);
h2 = (-B - sqrt(B.^2 - 4.*A.*C))./(2.*A);

f1 = matlabFunction(h1(s));
f2 = matlabFunction(h2(s));

x1t = @(t) real(f1(t));
y1t = @(t) imag(f1(t));
x2t = @(t) real(f2(t));
y2t = @(t) imag(f2(t));

figure
hold on
grid on
grid minor
axis equal
color = [0.8275 0.0745 0.7608];
plot(x1t(linspace(-pi, pi, 3000)), y1t(linspace(-pi, pi, 3000)), 'marker', '.', 'linestyle', 'none', 'color', color)
plot(x2t(linspace(-pi, pi, 3000)), y2t(linspace(-pi, pi, 3000)), 'marker', '.', 'linestyle', 'none', 'color', color)

%%


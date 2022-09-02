clear
% close all
clc

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

syms x w theta

%% find coefficients

% stencil = [0 1 2 3 4];

% stencil = 0;
% stencil = [0 1];
% stencil = [0 1 2];
% stencil = [0 1 3 4];
% stencil = [0 1 3 5 6];
% stencil = [0 1 3 5 7 8];
% stencil = [0 1 3 6 9 11 12];
stencil = [0 1 3 6 9 12 14 15];
% stencil = [0 1 3 6 9 12 15 17 18];
% stencil = [0 1 3 6 9 12 15 18 20 21];
% stencil = [0 1 3 6 9 13 17 20 23 25 26];
% stencil = [0 1 3 6 9 13 17 21 24 27 29 30];
% stencil = [0 1 3 6 10 14 18 22 26 29 31 32];
% stencil = [0 1 3 6 10 14 19 23 27 30 32 33];
% stencil = [0 1 3 6 9 13 17 21 25 28 31 33 34];
% stencil = [0 1 3 6 10 14 18 22 26 30 33 35 36];   % pretty
% stencil = [0 1 3 6 9 13 17 21 25 28 31 33 34];
% stencil = [0 1 3 6 10 14 19 23 28 32 36 39 41 42];
% stencil = [0 1 3 6 10 14 19 24 29 33 37 40 42 43];
% stencil = [0 1 3 6 10 14 19 24 29 34 38 42 45 47 48];
% stencil = [0 1 3 6 9 13 18 23 28 33 38 42 45 48 50 51]
% stencil = [0 1 3 6 10 14 19 24 29 34 39 43 47 50 52 53];
% stencil = [0 1 3 6 10 15 20 26 31 37 42 47 51 54 56 57];

n = length(stencil);
m = max(stencil);

stencil = sym(stencil);

lagrpol = @(j) prod((x - stencil(1:end ~= j)')./(stencil(j) - stencil(1:end ~= j)'), 1);
for i = 1:n
    if n == 1
        l = 1;
    else
        l(i) = lagrpol(n-i+1);
    end
end

alpha = [zeros(1, m), -1, 1];
betai = int(l, x, -1, 0);

for i = 1:n
    beta(stencil(i)+1) = betai(i);
end
beta = [beta 0];
% beta = [flip([1.99846268500208 -1.28551740667502 -0.217557223320753 0.797815926658251 -0.293203981664557]), 0]
k = length(beta)-1;

h(w) = sum(alpha.*w.^(0:k))./sum(beta.*w.^(0:k));

s(theta) = exp(1i*theta);

f = matlabFunction(h(s));

xt = @(t) real(f(t));
yt = @(t) imag(f(t));

% figure
hold on
grid on
grid minor
axis equal
color = lines(7);
plot(xt(linspace(-pi, pi, 3000)), yt(linspace(-pi, pi, 3000)))

% pretty(flip(beta(1:end-1)));

%%


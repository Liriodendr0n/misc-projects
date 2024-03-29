clear
close all
clc

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

syms x w theta

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
% stencil = [0 1 3 6 9 13 17 21 25 28 31 33 34];
% stencil = [0 1 3 6 10 14 18 22 26 30 33 35 36];
% stencil = [0 1 3 6 10 14 19 24 29 33 37 40 42 43];
% stencil = [0 1 3 6 10 14 19 24 29 34 38 42 45 47 48];
% stencil = [0 1 3 6 10 14 19 24 29 34 39 43 47 50 52 53];
% stencil = [0 1 2 3 4 5 6 7 8 9 10 11 12];

% stencil = [0 1 2 3 4];

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
betai = int(l, x, 0, 1);

for i = 1:length(stencil)
    beta(stencil(i)+1) = betai(i);
end
beta = [0 beta];
k = length(beta)-1;

h(w) = sum(alpha.*w.^(0:k))./sum(beta.*w.^(0:k));

s(theta) = exp(1i*theta);

f = matlabFunction(h(s));

xt = @(t) real(f(t));
yt = @(t) imag(f(t));

figure
hold on
grid on
grid minor
axis equal
color = lines(7);
plot(xt(linspace(-pi, pi, 3000)), yt(linspace(-pi, pi, 3000)), 'marker', '.', 'linestyle', 'none', 'color', color(2,:))
%%


clear
close all
clc

syms x w z(w) theta

%% find coefficients

% stencil = 0;
% stencil = [0 1];
% stencil = [0 1 2];
stencil = [0 1 2 3];
% stencil = [0 1 3 5 6];
% stencil = [0 1 4 7 10 11];
% stencil = [0 1 4 8 12 15 16];
% stencil = [0 1 3 6 9 12 15 18 20 21];
% stencil = [0 1 3 6 9 13 17 21 24 27 29 30];
% stencil = [0 1 3 6 10 14 18 22 26 30 33 35 36];

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
plot(xt(linspace(-pi, pi, 3000)), yt(linspace(-pi, pi, 3000)), 'marker', '.', 'linestyle', 'none', 'color', 'black')

%%


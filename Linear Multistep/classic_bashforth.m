clear
close all
clc

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

syms h x w z(w) phi xi yi t1 t2 t

%% find coefficients

stencil = [0 1 3 6 10 14 18 22 26 30 33 35 36];
%stencil = [0 1 2 3];
n = length(stencil);

lagrpol = @(j) prod((x - stencil(1:end ~= j)')./(stencil(j) - stencil(1:end ~= j)'), 1);
for i = 1:n
    if n == 1
        l = 2*x;
    else
        l(i) = lagrpol(n-i+1);
    end
end

bi = int(l, x, -1, 0);

for i = 1:length(stencil)
    b(stencil(i)+1) = bi(i);
end
n = length(b);
s(phi) = exp(1i*phi);
z(w) = (w.^n - w.^(n-1))/sum(b.*w.^(0:n-1));

f = matlabFunction(z(s));

xt = @(t) real(f(t));
yt = @(t) imag(f(t));

figure
hold on
grid on
grid minor
axis equal
fplot(xt, yt, [-pi pi])

xlim([-1.5, 0.5])
ylim([-1, 1])







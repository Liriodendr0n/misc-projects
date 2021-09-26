clear
close all
clc

syms h x w z(w) phi xi yi t1 t2 t

%% find coefficients

x0 = 1-[36 35 33 30 26 22 18 14 10 6 3 1 0];

n = length(x0);

lagrpol = @(j) prod((x - x0(1:end ~= j)')./(x0(j) - x0(1:end ~= j)'), 1);
for i = 1:n
    if n == 1
        l = 2*x;
    else
        l(i) = lagrpol(n-i+1);
    end
end

b = int(l, x, 0, 1);

s(phi) = exp(1i*phi);
z(w) = (w^n - w.^(n-1))/sum(b.*w.^(n:-1:1));

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







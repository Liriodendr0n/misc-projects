clear
close all
clc

%%

N = 100;

[xs, ys] = meshgrid(linspace(-1, 2, N), linspace(-1.25, 1.25, N));
xs = reshape(xs, N^2, 1);
ys = reshape(ys, N^2, 1);

kappa = 0;

f = @(t) kappa*t/2 .* (1-t);
fprime = @(t) kappa*(1/2 - t);
sigma = @(t) 1;
mu = @(t) 1;

PhiS = zeros(1, N^2);
PhiD = zeros(1, N^2);
for i = 1:N^2
    phiS = @(x, y) 1/(2*pi) * log(sqrt(x.^2 + y.^2));
    phiDx = @(x, y) 1/(2*pi) * -x./(x.^2 + y.^2);
    phiDy = @(x, y) 1/(2*pi) * -y./(x.^2 + y.^2);

    PhiS(i) = integral(@(t) sigma(t).*sqrt(1+fprime(t).^2) .* phiS(xs(i)-t, ys(i)-f(t)), 0, 1);
    PhiD(i) = integral(@(t) mu(t).*sqrt(1+fprime(t).^2) .* (cos(atan(fprime(t))).*phiDy(xs(i)-t, ys(i)-f(t)) - sin(atan(fprime(t))).*phiDx(xs(i)-t, ys(i)-f(t))), 0, 1);
end

PhiS = reshape(PhiS, N, N);
PhiD = reshape(PhiD, N, N);
xs = reshape(xs, N, N);
ys = reshape(ys, N, N);

figure
hold on
axis equal
view([0, 0, 1])
zlim([-1, 1])
contourf(xs, ys, PhiD)
fplot(f, [0, 1], 'color', 'black', 'linewidth', 2)
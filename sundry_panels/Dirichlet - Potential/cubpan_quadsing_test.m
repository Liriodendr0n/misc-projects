clear
close all
clc

%%

syms x y t real

beta1 = 1;
beta2 = 1;

phiDx(x, y) = -x/(x^2 + y^2);
phiDy(x, y) = -y/(x^2 + y^2);

z(t) = beta1/4 * (t+1)*(1-t)^2 + beta2/4 * (t+1)^2*(t-1);

zp = diff(z, t);

ds = simplify(sqrt(1 + diff(z, t, 1)^2));

mu(t) = sym(1) * t;

integrand = mu*ds * (cos(atan(zp))*phiDy(x-t, y-z(t)) - sin(atan(zp))*phiDx(x-t, y-z(t)));

fint = matlabFunction(integrand);

phiD = @(x, y) 1/(2*pi) * integral(@(t) fint(t, x, y), -1, 1, 'arrayvalued', false);

figure
fsurf(phiD, [-2, 2, -2, 2], 'edgecolor', 'none')
view(2)

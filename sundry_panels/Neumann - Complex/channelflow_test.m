clear
close all
clc

syms x y z t zeta

% zeta0 = (0);
% zeta1 = (0.7071+0.7071i);
syms zeta0 zeta1

nu = angle(zeta1-zeta0);

h = 1+1i;
% syms h

% f0 = 0;
% f1 = 1;
syms f0 f1

f(zeta) = f0 + (f1-f0)/(zeta1-zeta0) * (zeta - zeta0);

w(z, zeta) = 1/(2*h)*coth(pi/(2*h)*(z-zeta)) + 0/(2*h)*tanh(pi/(2*h)*(z-conj(zeta)));

W(z) = log(sinh(pi/(2*h)*(z-zeta1))/sinh(pi/(2*h)*(z-zeta0))) + log(cosh(pi/(2*h)*(z-conj(zeta1)))/cosh(pi/(2*h)*(z-conj(zeta0))));


w(z, zeta)*f(zeta)

figure
hold on
view(2)
fsurf(abs(w(x+1i*y, 0)), 'edgecolor', 'none', 'meshdensity', 100)
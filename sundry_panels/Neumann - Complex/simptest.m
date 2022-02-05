clear
close all
clc

%%

syms z zeta zeta0 zeta1 tau

f1 = -2i*tau*(1i*pi*zeta*log(1 + exp(i*pi*zeta/tau - i*pi*z/tau)) + tau*polylog(2, -exp(i*pi*zeta/tau - i*pi*z/tau))) + ...
     -2i*tau*(1i*pi*zeta*log(1 - exp(i*pi*zeta/tau - i*pi*z/tau)) + tau*polylog(2, exp(i*pi*zeta/tau - i*pi*z/tau)));

g1 = -2i*tau*(1i*pi*zeta*log(1 - exp(2i*pi/tau * (zeta - z))) + tau/2*polylog(2, exp(2i*pi/tau * (zeta - z))));

f2 = pi*(1i*zeta0*(1i*tau*log(1 + exp(1i*pi*z/tau - 1i*pi*zeta/tau)) + 1i*tau*log(1 - exp(1i*pi*z/tau - 1i*pi*zeta/tau))) + ...
    -tau*zeta0*log(exp(1i*pi*zeta/tau) + exp(1i*pi*z/tau)) - tau*zeta0*log(exp(1i*pi*zeta/tau) - exp(1i*pi*z/tau)) + ...
     pi*(2i*zeta0*zeta - 2i*zeta^2));
 
g2 = pi*(1i*zeta0*(1i*tau*log(1 - exp(2i*pi/tau * (z - zeta)))) - tau*zeta0*log(exp(2i*pi/tau * zeta) - exp(2i*pi/tau * z)) + 2i*pi*zeta*(zeta0 - zeta));


h = 4*pi^2*tau*(zeta0 - zeta1);

I(z, zeta, zeta0, zeta1) = (g1 + g2)/h;

Gsc(z, zeta0, zeta1) = I(z, zeta1, zeta0, zeta1) - I(z, zeta0, zeta0, zeta1); 

plotfun = matlabFunction(real(subs(Gsc(z, 0, 1), tau, 2+2i)));

[xs, ys] = meshgrid(linspace(-2, 2, 40), linspace(-2, 2, 40));

for i = 1:40
    i
    for j = 1:40
        zs(i, j) = plotfun(xs(i,j) + 1i*ys(i, j));
    end
end

figure
hold on
view(2)
surf(xs, ys, zs);

% figure
% hold on
% view(2)
% fsurf(@(x, y) abs(subs(Gsc(x+1i*y, 0, 1), tau, 2+2i)), 'meshdensity', 10, 'edgecolor', 'none')
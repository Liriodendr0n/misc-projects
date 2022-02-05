clear
close all
clc

%%

syms z z0 z1

% vortex panel
% w(z, z0, z1) = 1i/(2*pi)*log((z0-z)/(z1-z));

%% geometry

alpha = deg2rad(10);

N = 1000;


[x, y] = genNACA4([2 4 12], N+1);

ds = sqrt(diff(x).^2 + diff(y).^2);
nu = atan2(diff(y), diff(x))-pi/2;
xm = x(1:N) + diff(x)/2;
ym = y(1:N) + diff(y)/2;

s = [0; cumsum(ds)];

thetaTE = mod(nu(N)-nu(1)-pi, pi);

% figure
% hold on
% axis equal
% plot(x, y)
% quiver(xm, ym, cos(nu), sin(nu))
% % quiver(x, y, cos(nuj), sin(nuj))

%% influence matrices

z = xm + 1i*ym;
zeta = x + 1i*y;

A = imag(exp(1i*nu).*1i/(2*pi).*log((zeta(1:N).'-z)./(zeta(2:N+1).'-z)));
A(1:N+1:end) = -1/2;

A(N,:) = [1, zeros(1, N-2), 1];

b = real(exp(1i*nu).*exp(-1i*alpha));
b(N, 1) = 0;

g = A\b;

Cp = 1 - (real(exp(1i*(alpha+nu-pi/2)) + g/2)).^2;

figure
hold on
set(gca, 'ydir', 'reverse')
plot(xm, Cp)

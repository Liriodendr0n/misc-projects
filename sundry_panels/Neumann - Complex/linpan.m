clear
close all
clc

%% geometry

alpha = deg2rad(0);

N = 30;
wU = [0.16 0.16 0.22 0.12];
wL = [-0.09 -0.06 -0.12 0.12];
% set to zero
dzTE = 0;

[x, y] = genNACA4([2 4 12], N+1);
% [x, y] = genCST(wU, wL, dzTE, N+1);

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

% kth segment, jth endpoint

% off diagonal
K = [exp(-1i*nu.')/(2i*pi) .* (1 + (z-zeta(2:N+1).')./(zeta(2:N+1).'-zeta(1:N).') .* log((zeta(2:N+1).'-z)./(zeta(1:N).'-z))), zeros(N, 1)];
% diagonals
K(1:N+1:end) = (1/4 + 1/(2i*pi)) * exp(-1i*nu.');

% off diagonal
L = [zeros(N, 1), -exp(-1i*nu.')/(2i*pi) .* (1 + (z-zeta(1:N).')./(zeta(2:N+1).'-zeta(1:N).') .* log((zeta(2:N+1).'-z)./(zeta(1:N).'-z)))];
% diagonals
L(N+1:N+1:end) = (1/4 - 1/(2i*pi)) * exp(-1i*nu.');

C = K + L;

A(1, :) = [1, zeros(1, N)];
b(1, 1) = 0;

A(2:N+1, :) = imag(exp(1i*nu) .* C);
b(2:N+1, 1) = real(exp(1i*nu) .* exp(-1i*alpha));

A(N+2, 1:N+1) = [zeros(1, N), 1];
b(N+2, 1) = 0;

%% Solve

g = A\b;

Cp = 1 - g.^2;

Cl = 2*sum((g(1:N)+diff(g)/2).*ds);
disp("Cl = " + Cl)

% figure
% hold on
% set(gca, 'xdir', 'reverse')
% plot(s, g)

% figure
% hold on
% set(gca, 'ydir', 'reverse')
% plot(x, Cp)



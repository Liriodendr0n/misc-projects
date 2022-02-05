clear
close all
clc

%%

syms x y z zeta

% left and right endpoints
z0 = -1;
z1 = 1;

% real part = tangent doublet strength (source-like)
% imag part = normal doublet strength (vortex-like)
f0 = 1i;
f1 = 1i;

% complex potentials (constant and linear strength)
phi0(z) = 1/sym(2*pi) * f0 * log((z-z1)/(z-z0));
phi1(z) = 1/sym(2*pi) * (f0 + ((f0-f1)*(z-z0))/(z0-z1)) * log((z-z1)/(z-z0));


figure
hold on
view(2)
zlim([-5, 5])
fsurf(real(phi1(x+1i*y)), 'edgecolor', 'none')



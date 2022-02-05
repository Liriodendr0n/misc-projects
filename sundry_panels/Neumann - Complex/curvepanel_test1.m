clear
close all
clc

%%

syms z zeta

% path
gamma(zeta) = zeta + 1i*zeta^3 + 1i*zeta^2 + 1i*zeta;

% strength
f(zeta) = 1+zeta+zeta^2;

% green's function
% G(z, zeta) = 1/sym(2*pi) * 1/(z - zeta);
G(z, zeta) = 1/sym(2*pi) * log(z - zeta);

% integrand
G(z, f(zeta)) * diff(gamma, zeta) * f(zeta)

% I(z, zeta) = -(1i*log(zeta^2-1i*zeta+1i*z)+(4*1i*z+1)*sqrt(-1/(4*1i*z+1))*log(2*zeta+sqrt(-4*1i*z-1)-1i)+(-4*1i*z-1)*sqrt(-1/(4*1i*z+1))*log(2*zeta-sqrt(-4*1i*z-1)-1i)+4*zeta)/(4*pi);


% zeta0 = -1;
% zeta1 = 1;
% 
% IF(z) = simplify(I(z, zeta1) - I(z, zeta0));
% 
% syms x y
% figure
% hold on
% view(2)
% fsurf(imag(IF(x+1i*y)), 'edgecolor', 'none')
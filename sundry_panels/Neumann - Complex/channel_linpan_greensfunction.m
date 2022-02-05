clear
close all
clc

%% linear strength (increasing from 0 to 1 from z0 to z1)

syms z zeta

f0 = 0;
f1 = 1;

syms z0 z1
% z0 = 0; % must be 0 for the antiderivative
% z1 = 1+0.5i;

nu = angle(z1-z0);

f(zeta) = f0 + (f1-f0)/(z1-z0) * (zeta - z0);
Z(zeta) = z0 + (z1-z0)*zeta;

% c = sym(pi)/4;
c = 1;

w(z, zeta) = coth(c*(z - zeta));

% I(z, zeta) = (polylog(2, 1-exp(2*c*(zeta-z))) - 2*c*z*log(exp(2*c*zeta)-exp(2*c*z)) + c^2*zeta^2)/(2*c^2*(z1-z0));
% I(z, zeta) = (polylog(2, 1-exp(2*c*(zeta-z))) - 2*c*z*log(exp(2*c*zeta)-exp(2*c*z)) + c^2*zeta^2)/(2*c^2*(z1-z0));

W(z) = (2*c*z0*log( sinh(c*(z1-z)) / (sinh(c*(z0-z))) ) + polylog(2, 1-exp(2*c*(z1-z))) - polylog(2, 1-exp(2*c*(z0-z))) + 2*c*z*log((exp(2*c*z0)-exp(2*c*z))/(exp(2*c*z1)-exp(2*c*z))) + c^2*z1^2 - c^2*z0^2)/(2*c^2*(z1-z0));

% W2(z) = log((exp(2*c*z0)-exp(2*c*z))/(exp(2*c*z1)-exp(2*c*z))) + log( sinh(c*(z1-z)) / (sinh(c*(z0-z))) ) + c;
% W3(z) = log( sinh(c*(z1-z)) / (sinh(c*(z0-z))) );

% W3(z) = dilog((z-z0)*exp(-1i*(nu-pi))) - dilog((z1-z)*exp(-1i*(nu)));

%%

% syms x y
% 
% figure
% hold on
% view(2)
% fsurf(imag(W3(x+1i*y)), 'edgecolor', 'none', 'meshdensity', 10)


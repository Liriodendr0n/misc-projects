clear
close all
clc


%% simpler form (katz & plotkin)

syms x y theta1 theta2 r1 r2 real

x1 = 0;
x2 = 1;
theta1 = atan((x-x1)/y);
theta2 = atan((x-x2)/y);
r1 = (x-x1)^2 + y^2;
r2 = (x-x2)^2 + y^2;

mu0 = 1;
mu1 = 0;
mu2 = 0;

sigma0 = 0;
sigma1 = 0;

Phi012(x, y) = - mu0/(2*sym(pi))*(theta1 - theta2) - mu1/(4*sym(pi))*(2*x*(theta1 - theta2) + y*log(r2/r1)) + mu2/(2*sym(pi)) * ((x^2-y^2)*(theta2-theta1) - x*y*log(r2/r1) + y*(x1-x2));
Phi012_p(x) = sym(0);
% doublet potential is zero on y = 0
% potential jump across y = 0 is mu0 + mu1x + mu2x^2

Sigma01(x, y) = sigma0/(4*sym(pi)) * ((x-x1)*log(r1) - (x-x2)*log(r2) + 2*y*(theta1 - theta2)) + sigma1/(8*sym(pi)) * ((x^2-x1^2-y^2)*log(r1) - (x^2-x2^2-y^2)*log(r2) + 4*x*y*(theta1-theta2) - 2*x*(x2-x1));
Sigma01_p(x) = sigma0/(4*sym(pi)) * (x*log(x^2) - log((x - 1)^2)*(x - 1)) - sigma1/(8*sym(pi)) * (2*x + log((x - 1)^2)*(x^2 - 1) - x^2*log(x^2));

% source potential is zero at linear endpoints
% 0 at 0 end of linear, -1/(4pi) at 1 end of linear


%% Field Potential

theta = 30;

figure
hold on
axis equal
fsurf(Phi012+Sigma01, [-1 2 -1.5 1.5], 'Edgecolor', 'none')
% fcontour(atan2(y,(1-x)/(2*pi)), 'fill', 'on')
view([0, 0, 1])
xlim([-1 2])
ylim([-1.5 1.5])
zlim([-1.5 1.5])
caxis([-0.5, 0.5])

%% Panel Potential

% figure
% hold on
% grid on
% fplot(Sigma01_p, [-1, 2])

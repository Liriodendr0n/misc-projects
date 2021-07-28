clear
close all
clc

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% chamber properties

g = 7/5;
T0 = 3000;
R = 8314.5/2;

a0 = sqrt(g*R*T0);

p0 = 1;
rho0 = p0/(R*T0);

Nchar = 100;

%% speeds of sound

a = linspace(a0, 0, Nchar);

u = 2/(g-1) * (a0 - a);

T = T0*(a./a0).^2;
p = p0*(T./T0).^(g/(g-1));
rho = rho0*(T./T0).^(1/(g-1));

%% characteristic slopes

dtdx = 1./(u-a);


%% end of simple region
x0 = 1;
xmax = x0*(2/(g-1) + 1);
tEnd = -x0*dtdx(1);


%% wavelet locations

xw = (tEnd + dtdx*x0)./dtdx;


%% plot test

figure('name', 'thermoprop', 'units', 'inches', 'papersize', [6 6])
hold on
grid on
grid minor
% for plotind = 1:Nchar
%     fplot(@(x) dtdx.*(x - x0), [0, xmax])
% end
% 
xlim([0, xmax])
%ylim([0, tEnd])
% scatter(xw, ones(size(xw))*tEnd)
%plot(xw, u)
plot(xw, u./u(end))
plot(xw, a./a0)
plot(xw, T./T0)
plot(xw, rho./rho0)
plot(xw, p./p0)
plot(xw, (rho.*u./(rho0*u(end)))./max(rho.*u./(rho0*u(end))), 'linestyle', '--')
plot(xw, (rho.*u.^2./(rho0*u(end).^2))./max(rho.*u.^2./(rho0*u(end).^2)), 'linestyle', '--')
legend('$u/u_{max}$', '$a/a_0$', '$T/T_0$', '$\rho/\rho_0$', '$p/p_0$', '$(\rho u)/(\rho u)_{max}$', '$(\rho u^2)/(\rho u^2)_{max}$', 'location', 'east')
title('Simple Region Thermodynamic Properties')
xlabel('Length [m]')
ylabel('Relative Value')

print('SimpleTube.png', '-dpng', '-r300')
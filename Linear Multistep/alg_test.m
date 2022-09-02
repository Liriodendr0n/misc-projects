clear
close all
clc

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% model problem
% exponential for checking stability
lambda = 1i;
f = @(t, y) lambda*y;
y0 = 1;

% % polynomial for checking order (highest exact solution is order)
% m = 7;
% f = @(t, y) m*t.^(m-1); 
% y0 = 1;

% % % r theta orbit IC = (rdot thetadot r theta)
% f = @(t, y) [y(3).*y(2).^2 - 1./y(3).^2; -2./y(3) .*y(1).*y(2); y(1); y(2)];
% y0 = [0; 1.1; 1; 0];

% % x y orbit IC = (vx vy rx ry)
% f = @(t, y) [-y(3)./((y(3).^2 + y(4).^2).^(3/2)); -y(4)./((y(3).^2 + y(4).^2).^(3/2)); y(1); y(2)];
% y0 = [0; 0.5; 1; 0];

tspan = [0 10];

h = 0.1;

%% stencil and start

% stencil = [0 1 3 4];                                        % Usable 4th order
stencil = [0 1 3 6 9 12 14 15];                             % Usable 8th order
% stencil = [0 1 3 6 10 14 19 23 27 30 32 33];                % Usable 12th order
% stencil = [0 1 3 6 10 15 20 26 31 37 42 47 51 54 56 57];    % Usable 16th order
% stencil = [0 1 3 6 10 15 21 27 33 40 47 54 60 66 72 77 81 84 86 87];

%% test integrator

% tic
% [t, y] = bashforth(f, tspan, h, stencil, y0);
% toc


[t, y] = explicit_leastSquares(f, tspan, h, 20, 10, y0);


%% plots

% for checking stability
figure
hold on
grid on
grid minor
fplot(@(t) real(exp(lambda*t)*y0), tspan)
plot(t, real(y), 'marker', '.')
title("$h\lambda$ = " + h*lambda)

yE = exp(lambda*(t(end)-t(1)))*y0;
yA = y(end);
errMag = log10(abs((yA - yE)/yA))


% % for checking order with f = @(t, y) m*t.^(m-1);
% figure
% hold on
% grid on
% grid minor
% set(gca, 'Yscale', 'log')
% plot(t, abs(y - t.^(m))-1)
% ylim([1e-17 1e17])

% % for plotting orbits
% figure
% hold on
% axis equal
% grid on
% grid minor
% % plot(y(3,:).*cos(y(4,:)), y(3,:).*sin(y(4,:)))
% plot(y(3,:), y(4,:))%, 'marker', '.')
% % xlim([-2 2])
% % ylim([-2 2])
% % % orbital energy xy
% E0 = -1./sqrt(y0(3).^2+y0(4).^2) + 1/2*(y0(1).^2+y0(2).^2);
% E = -1./sqrt(y(3,:).^2+y(4,:).^2) + 1/2*(y(1,:).^2+y(2,:).^2);
% figure
% plot(t, log10(abs(E-E0)))


clear
close all
clc

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% model problem
% polynomial for checking order
% f = @(t, y) m*t.^(m-1); 

% r theta orbit IC = (rdot thetadot r theta)
% f = @(t, y) [y(3).*y(2).^2 - 1./y(3).^2; -2./y(3) .*y(1).*y(2); y(1); y(2)];

% x y orbit IC = (vx vy rx ry)
f = @(t, y) [-y(3)./((y(3).^2 + y(4).^2).^(3/2)); -y(4)./((y(3).^2 + y(4).^2).^(3/2)); y(1); y(2)];

y0 = [0; 1.2; 1; 0];
tspan = [0 15];

h = 0.01;

%% stencil and start

stencil = [0 1 3 6 10 14 18 22 26 30 33 35 36];
%stencil = [0 1 3 5 6];

tstart = -(0:max(stencil))*h;

opts = odeset('reltol', 1e-12);
[t2, y2] = ode45(@(t, y) f(-t, y), -tstart, y0);
start45 = ode45(@(t, y) f(-t, y), [0 (max(stencil)*h)], y0, opts);
ystart = deval(start45, -tstart);


%% test integrator
tic
[t, y] = bashforth(f, tspan, h, stencil, ystart);
toc
%% plots

% % for checking order with f = @(t, y) m*t.^(m-1);
% figure
% hold on
% grid on
% grid minor
% set(gca, 'Yscale', 'log')
% plot(t, abs(y - t.^(m)))
% ylim([1e-17 1e17])

figure
hold on
axis equal
grid on
grid minor
% plot(y(3,:).*cos(y(4,:)), y(3,:).*sin(y(4,:)))
plot(y(3,:), y(4,:))
xlim([-2 2])
ylim([-2 2])


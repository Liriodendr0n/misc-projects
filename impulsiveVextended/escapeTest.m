clear
close all
clc

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

tstop = 1e3;
tspan = [0, tstop];

% minimum acceleration for a radial out burn
%alpha = 0.09007521167227;
alpha = 0.126;

opts = odeset('RelTol',1e-12, 'AbsTol',1e-12, 'Events', @escapeV);
odeSoln = ode45(@(t, y) odefun(t, y, alpha), tspan, [0; 1; 1; 0], opts);

%tesc = tstop;
tesc = odeSoln.xe

t = linspace(0, tesc, 1000);
y = deval(odeSoln, t);

t = t;

rdot = y(1,:);
thetadot = y(2,:);
r = y(3,:);
theta = y(4,:);

% r = odeSoln.y(3,:);
% theta = odeSoln.y(4,:);

figure('name', 'spiral', 'units', 'inches', 'papersize', [6 6], 'paperposition', [0 0 6 6])
hold on
grid on
plot(r.*cos(theta), r.*sin(theta))
plot(cos(linspace(0, 2*pi, 360)), sin(linspace(0, 2*pi, 360)))
axis equal

xlim([floor(min(r.*cos(theta))), ceil(max(r.*cos(theta)))])
ylim([floor(min(r.*sin(theta))), ceil(max(r.*sin(theta)))])

xlabel('x/R')
ylabel('y/R')


%print('spiral.pdf', '-dpdf', '-painters')

% figure
% plot(t, atan2d(rdot, r.*thetadot))

Dv = alpha*tesc

T = 0.5 * (rdot.^2 + (r.*thetadot).^2);
V = -1./r;

E = T + V;
L = T - V;

v = sqrt(rdot.^2 + (r.*thetadot).^2);

% figure
% plot(t, r)

figure
plot(t, v./sqrt(1./r))


%% orbit integrator function, polar coordinates
function dydt = odefun(t, y, alpha)
    
dydt = zeros(4, 1);
    
% d/dt (rdot) = r * thetadot^2 - mu/r^2 + alpha*sin(arctan(rdot/thetadot))
%dydt(1) = y(3)*y(2)^2 - 1/y(3)^2 + alpha*sin(atan2(y(1),y(2)*y(3)));
dydt(1) = y(3)*y(2)^2 - 1/y(3)^2 + alpha;

% d/dt thetadot = -2/r * rdot * thetadot + alpha*cos(arctan(rdot/thetadot))
%dydt(2) = -2/y(3) *y(1)*y(2) + 1/y(3) *alpha*cos(atan2(y(1),y(2)*y(3)));
dydt(2) = -2/y(3) *y(1)*y(2);
    
% d/dt r = rdot
dydt(3) = y(1);

% d/dt theta = thetadot
dydt(4) = y(2);

% get rid of the stupid highlight
t = t;

end

%% stop when escape velocity reached
function [value,isterminal,direction] = escapeV(t, y)

value = sqrt(2/y(3)) - sqrt(y(1)^2 + (y(2)*y(3))^2);
isterminal = 1;
direction = 0;

end

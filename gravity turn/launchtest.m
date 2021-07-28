clear
close all
clc


%% parameters and initial conditions

% mu = 1;

alpha = 2;

v0 = 1e-6;
gamma0 = fsolve(@(gamma0) launchODE(alpha, v0, gamma0), 1);

r0 = 1;
theta0 = 0;
rdot0 = v0*tan(gamma0)*sqrt(1/(tan(gamma0)^2 +1));
thetadot0 = v0/r0*sqrt(1/(tan(gamma0)^2 +1));

y0 = [rdot0; thetadot0; r0; theta0];

tstop = 1e3;
tspan = [0, tstop];

opts = odeset('RelTol',1e-12, 'AbsTol',1e-12, 'Events', @escapeV);
odesoln = ode45(@(t, y) odefun(t, y, alpha), tspan, y0, opts);

tend = odesoln.xe;
yend = odesoln.ye;

rdot = odesoln.y(1,:);
thetadot = odesoln.y(2,:);
r = odesoln.y(3,:);
theta = odesoln.y(4,:);

t = odesoln.x;
v = sqrt(rdot.^2 + (r.*thetadot).^2);
gamma = atan2(rdot, r.*thetadot);
rddot = r.*thetadot.^2 - 1./(r.^2) + alpha.*sin(gamma);
thetaddot = -2./r .*rdot.*thetadot + alpha./r .*cos(gamma);

rdotend = yend(1);
thetadotend = yend(2);
rend = yend(3);
thetaend = yend(4);

vend = sqrt(rdotend^2 + (rend*thetadotend)^2);
gammaend = atan2(rdotend, rend*thetadotend)

DvLaunch = alpha*tend
rend-1

%% plots

% trajectory
figure
hold on
axis equal
plot(ones(1,360*4).*cos(linspace(0, 2*pi, 360*4)+pi/2), ones(1,360*4).*sin(linspace(0, 2*pi, 360*4)+pi/2))
plot(odesoln.y(3,:).*cos(odesoln.y(4,:)+pi/2), odesoln.y(3,:).*sin(odesoln.y(4,:)+pi/2))
%plot(rend*ones(1,360*4).*cos(linspace(0, 2*pi, 360*4)+pi/2), rend*ones(1,360*4).*sin(linspace(0, 2*pi, 360*4)+pi/2))


% launch data
figure
hold on
grid on
plot(t, r)
plot(t, rdot)
plot(t, rddot)

figure
hold on
axis equal
plot(theta, r-1)

%% orbit integrator function, polar coordinates
function dydt = odefun(t, y, alpha)
    
dydt = zeros(4, 1);
    
% d/dt (rdot) = r * thetadot^2 - mu/r^2 + alpha*sin(arctan(rdot/thetadot))
dydt(1) = y(3)*y(2)^2 - 1/y(3)^2 + alpha*sin(atan2(y(1),y(2)*y(3)));

% d/dt thetadot = -2/r * rdot * thetadot + alpha*cos(arctan(rdot/thetadot))
dydt(2) = -2/y(3) *y(1)*y(2) + 1/y(3) *alpha*cos(atan2(y(1),y(2)*y(3)));

% d/dt r = rdot
dydt(3) = y(1);

% d/dt theta = thetadot
dydt(4) = y(2);

% get rid of the stupid highlight
t = t;

end

%% stop when escape velocity reached
function [value,isterminal,direction] = escapeV(t, y)

value = sqrt(1/y(3)) - sqrt(y(1)^2 + (y(2)*y(3))^2);
isterminal = 1;
direction = 0;

end
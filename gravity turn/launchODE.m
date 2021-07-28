function [gammaend] = launchODE(alpha, v0, gamma0)
%LaunchODE integrates equations of motion during a rocket launch
%   Detailed explanation goes here

%% setup and results
r0 = 1;
theta0 = 0;
rdot0 = v0*tan(gamma0)*sqrt(1/(tan(gamma0)^2 +1));
thetadot0 = v0/r0*sqrt(1/(tan(gamma0)^2 +1));

y0 = [rdot0; thetadot0; r0; theta0];

tstop = 10;
tspan = [0, tstop];

opts = odeset('RelTol',1e-12, 'AbsTol',1e-12, 'Events', @escapeV);
odesoln = ode45(@(t, y) odefun(t, y, alpha), tspan, y0, opts);

tend = odesoln.xe;
yend = odesoln.ye;

rdotend = yend(1);
thetadotend = yend(2);
rend = yend(3);
%thetaend = yend(4);

%vend = sqrt(rdotend^2 + (rend*thetadotend)^2);
gammaend = atan2(rdotend, rend*thetadotend);

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

end

%% stop when orbital velocity reached
function [value,isterminal,direction] = escapeV(t, y)

value = sqrt(1/y(3)) - sqrt(y(1)^2 + (y(2)*y(3))^2);
isterminal = 1;
direction = 0;

end

end


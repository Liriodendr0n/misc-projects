function [Dv, tesc, gammaesc, thetaesc] = escDv_tan(alpha)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


tstop = 1e5;
tspan = [0, tstop];


opts = odeset('RelTol',1e-6, 'AbsTol',1e-6, 'Events', @escapeV);
odeSoln = ode45(@(t, y) odefun(t, y, alpha), tspan, [0; 1; 1; 0], opts);

tesc = odeSoln.xe;

% t = linspace(0, tesc, tstop);
% y = deval(odeSoln, t);

rdot = odeSoln.y(1,:);
thetadot = odeSoln.y(2,:);
r = odeSoln.y(3,:);
theta = odeSoln.y(4,:);

% r = odeSoln.y(3,:);
% theta = odeSoln.y(4,:);

% figure
% plot(r.*cos(theta), r.*sin(theta))
% axis equal
% 
% figure
% plot(t, atan2d(rdot, r.*thetadot))

Dv = alpha*tesc;
gammaesc = atan2(rdot(end), r(end)*thetadot(end));

thetaesc = theta(end);

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

value = sqrt(2/y(3)) - sqrt(y(1)^2 + (y(2)*y(3))^2);
isterminal = 1;
direction = 0;

end



end


clear
close all
clc


tstop = 10000;
tspan = [0, tstop];

opts = odeset('RelTol',1e-6, 'AbsTol',1e-6);
odeSoln = ode45(@odefun, tspan, [0; 1.4; 1; 0], opts);

t = linspace(0, tstop, tstop*10);
y = deval(odeSoln, t);

rdot = y(1,:);
thetadot = y(2,:);
r = y(3,:);
theta = y(4,:);

% r = odeSoln.y(3,:);
% theta = odeSoln.y(4,:);


% plot(r.*cos(theta), r.*sin(theta))
% axis equal

plot(t, thetadot)



%% orbit integrator function, polar coordinates
function dydt = odefun(t, y)
    
dydt = zeros(4, 1);
    
% d/dt (rdot) = r * thetadot^2 - mu/r^2
dydt(1) = y(3)*y(2)^2 - 1/y(3)^2;
    
% d/dt thetadot = -2/r * rdot * thetadot
dydt(2) = -2/y(3) *y(1)*y(2);
    
% d/dt r = rdot
dydt(3) = y(1);

% d/dt theta = thetadot
dydt(4) = y(2);

% get rid of the stupid highlight
t = t;

end
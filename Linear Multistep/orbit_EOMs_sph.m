clear
close all
clc

syms r rdot rddot theta thetadot thetaddot phi phidot phiddot real
syms mu J2 real

% state vector
x = [r; theta; phi; rdot; thetadot; phidot];

%% inertial coordinates
% [r radius; theta longitude; phi colatitude]
% R, V, A from illinois dynref (Matthew West)
R = [r; 0; 0];
V = [rdot; r*thetadot*sin(phi); r*phidot];
A = [rddot - r*thetadot^2*sin(phi)^2 - r*phidot^2; r*thetaddot*sin(phi) + 2*rdot*thetadot*sin(phi) + 2*r*thetadot*phidot*cos(phi); r*phiddot + 2*rdot*phidot - r*thetadot^2*sin(phi)*cos(phi)];

eqn = A == [-mu/(r^2); 0; 0];

soln = solve(eqn, [rddot, thetaddot, phiddot]);
EOMs(x) = [rdot; thetadot; phidot; soln.rddot; soln.thetaddot; soln.phiddot];

%% rotation about z (phi = 0)
syms omega real
Omega = [omega*cos(phi); 0; -omega*sin(phi)];

Vr = V - cross(Omega, R);
Ar = A - cross(2*Omega, Vr) - cross(Omega, cross(Omega, R));

eqnr = Ar == [-mu/(r^2); 0; 0];

solnr = solve(eqnr, [rddot, thetaddot, phiddot]);
EOMsr(x) = [rdot; thetadot; phidot; solnr.rddot; solnr.thetaddot; solnr.phiddot];

%% linearization

M(x) = jacobian([solnr.rddot; solnr.thetaddot; solnr.phiddot], [rdot; thetadot; phidot]);

%%

w = 0;

numEOMs = matlabFunction(subs(EOMsr, [mu omega], [1 w]));

r0 = 1;
theta0 = 0;
phi0 = pi/2;
rdot0 = 0.2;
thetadot0 = 1;
phidot0 = 0.4;

tspan = [0 2*pi]*16;
teval = linspace(tspan(1), tspan(2), 1000);
opts = odeset('reltol', 1e-6);
odesoln = ode45(@(t, x) numEOMs(x(1), x(2), x(3), x(4), x(5), x(6)), tspan, [r0; theta0; phi0; rdot0; thetadot0; phidot0], opts);
odevals = deval(odesoln, teval);

ts = teval;
rs = odevals(1,:);
thetas = odevals(2,:);
phis = odevals(3,:);
rdots = odevals(4,:);
thetadots = odevals(5,:);
phidots = odevals(6,:);


figure
hold on
grid on
grid minor
axis equal
xlim([-1.5 1.5])
ylim([-1.5 1.5])
zlim([-1.5 1.5])
view(3)
plot3(rs.*cos(thetas).*sin(phis), rs.*sin(thetas).*sin(phis), rs.*cos(phis))
scatter3(0, 0, 0, 1000, 'marker', '.')

% figure
% hold on
% grid on
% grid minor
% 
% plot(ts, [rdots; thetadots; phidots])

%%


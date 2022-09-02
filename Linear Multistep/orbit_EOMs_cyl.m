clear
close all
clc

syms r rdot rddot theta thetadot thetaddot real
syms mu real
% r radius
% theta longitude

R = [r; 0];
V = [rdot; r*thetadot];
A = [rddot - r*thetadot^2; r*thetaddot + 2*rdot*thetadot];

eqn = A == [-mu/(r^2); 0];

soln = solve(eqn, [rddot, thetaddot]);
EOMs = [soln.rddot; soln.thetaddot];


clear
close all
clc

%%

syms r theta rdot thetadot real

x = [r, theta, rdot, thetadot];
f = [rdot; thetadot; r*thetadot^2 - 1/r^2; -2/r * rdot * thetadot];

J = jacobian(f, x)
f
x;

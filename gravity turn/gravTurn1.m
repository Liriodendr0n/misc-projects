clear
close all
clc


%%

alpha = 1.4;
v0 = 0.01;

gamma0 = rad2deg(fsolve(@(gamma0) launchODE(alpha, v0, gamma0), 1))
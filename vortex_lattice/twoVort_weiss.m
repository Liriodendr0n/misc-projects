clear
close all
clc

%%
syms rho u G1 G2 L M Pi

eqn1 = L/(rho*u) == (G1 + G2);
eqn2 = M/(rho*u) == G1*(sin(Pi/10)^2 - 1/4) + G2*(sin(3*Pi/10) - 1/4);

soln = solve([eqn1, eqn2], [G1, G2])

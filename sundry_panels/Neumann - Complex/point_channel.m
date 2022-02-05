clear
close all
clc

%%

syms z

h = 1i;
zeta = 0;

V(z) = cot(pi/h * (z-zeta))/(2*h) - 0/(2*pi*(z-zeta));

syms x y

figure
hold on
view(2)
fsurf(real(V(x+1i*y)), 'edgecolor', 'none')
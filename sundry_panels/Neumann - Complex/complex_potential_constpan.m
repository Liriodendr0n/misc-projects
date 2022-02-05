clear
close all
clc

%%

syms x y z

z0 = 0;
z1 = 1i;

% phi(z) = log(((z-z1)^z1)/((z-z0)^z0) * ((z-z0)/(z-z1))^z);

phi(z) = log(sinh(z));


figure
hold on
view(2)
fsurf(imag(phi(x+1i*y)))
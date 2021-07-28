clear
close all
clc

r_p = 0.4572/2;
A_p = pi*r_p^2;

m_p = 70;

p0 = 500*10^6;
T0 = 3000;
L_b = 10;

M = 0.002;
R = 8.3145/M;
g = 1.4;

a0 = sqrt(g*R*T0);

pfun = @(u) p0*(1 - (g-1)/2 * u/a0)^((2*g)/(g-1));

Afun = @(u) pfun(u)*A_p/m_p;

[t, v] = ode45(@(t,v) Afun(v), [0 0.01], 0);

r = cumtrapz(t, v);
a = gradient(v, t);

figure
plot(r, v)
xlim([0, L_b])

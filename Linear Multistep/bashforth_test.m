clear
close all
clc

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% model problem

m = 14;

f = @(t, y) m*t.^(m-1);

y0 = zeros(1, 1);
tspan = [0 3];

h = 0.03;

%% stencil and start

stencil = [0 1 3 6 10 14 18 22 26 30 33 35 36];
%stencil = [0 2 5 7];

tstart = -(0:max(stencil))*h;

ystart = exp(m.*tstart).*y0;

%% test integrator

[t, y, f] = bashforth(f, tspan, h, stencil, ystart);

%% plots

figure
hold on
grid on
grid minor
set(gca, 'Yscale', 'log')
plot(t, abs(y - t.^(m)))
ylim([1e-17 1e17])



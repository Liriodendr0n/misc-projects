clear
close all
clc

%% parameters

a = 1;
e = 0.5;

Theta = 1;

N = 101;
thetas = linspace(-pi, pi, N);


%% shape calculations

% angular function
phi = @(theta) 2*atan(sqrt((1+e)/(1-e)).*tan(invKep(theta, e)/2));
% phi = @(theta) 2*atan(sqrt((1+e)/(1-e)).*tan(theta/2));
phis = phi(thetas');

% transfer function
n = gradient(phis, thetas);

% gear shape
% radius 1
R1 = a*n./(1+n);

% radius 2
R2 = a./(1+n);

%% arc length reparameterization

% gear 1
ds1 = sqrt(gradient(R1, thetas).^2 + R1.^2);
s1 = cumtrapz(thetas, ds1);

s1invFun = @(theta) interp1(s1, thetas, theta, 'pchip');
s1inv = s1invFun(max(s1)*(thetas+pi)/(2*pi));

R1invFun = @(theta) interp1(thetas, R1, theta, 'pchip');
R1inv = R1invFun(s1inv);

% gear 2
ds2 = sqrt(gradient(R2, phi(thetas)).^2 + R2.^2);
s2 = cumtrapz(phi(thetas), ds2);

s2invFun = @(theta) interp1(s2, thetas, theta, 'pchip');
s2inv = s1invFun(max(s2)*(thetas+pi)/(2*pi));

R2invFun = @(theta) interp1(thetas, R2, theta, 'pchip');
R2inv = R2invFun(s2inv);

%% cartesian
x1s = R1inv.*cos(s1inv);
y1s = R1inv.*sin(s1inv);

x2s = R2inv.*cos(phi(s2inv));
y2s = R2inv.*sin(phi(s2inv));

s1s = std(sqrt(diff(x1s).^2 + diff(y1s).^2))
s2s = std(sqrt(diff(x2s).^2 + diff(y2s).^2))

%% plots

figure
hold on
grid on
grid minor
axis equal
plot(R1inv.*cos(s1inv+Theta), R1inv.*sin(s1inv+Theta), 'marker', '.')
plot(a - R2inv.*cos(phi(s2inv)+phi(Theta)), R2inv.*sin(phi(s2inv)+phi(Theta)), 'marker', '.')

plot([0 a], [0 0], 'linestyle', 'none', 'marker', '.', 'color', 'black')

xlim([-1, 2])
ylim([-1, 1])

%%


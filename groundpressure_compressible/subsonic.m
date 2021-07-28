clear
close all
clc

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% input parameters


L = 1;
Uinf = 1e+3;
H = 1000;

kz = L/Uinf;


% Nx = 1001;
% Ny = 1001;
% Nz = 3;

NR = 101;
NTheta = 101;
Nz = 1;

Nx = NR;
Ny = NR;

M = 0.0;

Beta = sqrt(1 - M^2);

%% perturbation velocity potential function
% doublet line form, far field

syms x y z beta h

phi_doublet = 1/(4*pi) * z./(y.^2+z.^2).*(1+x./(sqrt(x.^2+beta.^2.*(y.^2+z.^2))));

phi = kz.*subs(phi_doublet, z, z-h) + -kz.*subs(phi_doublet, z, z+h);
phiFun = matlabFunction(phi);

%% grid points

Rs = logspace(-3, 2, NR)*H;
Thetas = linspace(-pi, pi, NTheta)+1e-6;

xs = (Rs).*cos(Thetas').*ones(1, 1, Nz);
ys = (Rs./Beta).*sin(Thetas').*ones(1, 1, Nz);

zs = zeros(NTheta, NR, Nz)+0*H;

% scatter3(reshape(xs, [prod(size(xs)), 1]), reshape(ys, [prod(size(ys)), 1]), reshape(zs, [prod(size(zs)), 1]))

%% perturbation velocity at each point

u = diff(phi, x);
v = diff(phi, x);
w = diff(phi, x);

ufun = matlabFunction(u);
vfun = matlabFunction(v);
wfun = matlabFunction(w);

us = ufun(Beta, H, xs, ys, zs);
vs = vfun(Beta, H, xs, ys, zs);
ws = wfun(Beta, H, xs, ys, zs);

%% pressure perturbation at each point
Pp = -Uinf*us - us.^2 - vs.^2 - ws.^2;

%% check integrated lift

%Lcheck = trapz(Rs, Rs.*trapz(Thetas, Pp(:,:,1)))/beta;

%% make plots

figure('units', 'inches', 'papersize', [4 4])
hold on
axis equal
contourf(xs(:, :, 1)./H, ys(:, :, 1)./H, log10(abs(Pp(:, :, 1)./(L/(2*pi*H^2)))+1e-3), -3:0.01:1, 'edgecolor', 'none')
%contourf(reshape(xp(:, :, 2)./h, [Nx, Ny]), reshape(yp(:, :, 2)./h, [Nx, Ny]), reshape(Pp(:, :, 2), [Nx, Ny])./(L/(2*pi*h^2)), linspace(-1, 1, 30), 'edgecolor', 'none')
colorbar
c = colorbar;
c.Label.String = 'Log Pressure Disturbance';
c.Label.Interpreter = 'latex';
set(c,'TickLabelInterpreter','latex')
%colormap(parula(80))
caxis([-3 1])
xticks([-50:10:50])
yticks([-50:10:50])
xlabel('x/h')
ylabel('y/h')
xlim([-50, 50])
ylim([-50, 50])
title("Ground Pressure: M = "+ floor(M)+"."+floor(10*(M-floor(M)))+floor(10*(M*10-floor(M*10))))


%print("Mach"+M+".png", '-dpng', '-r600')

% figure
% hold on
% plot(Rs, Pp(1,:,2))
% plot(Rs, L/(2*pi*H^2).*cos(atan(Rs./H)).^3)




%%


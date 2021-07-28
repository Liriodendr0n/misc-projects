clear
close all
clc

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% input parameters


U = 1e3;
G = 1;
b = 1e-3;

L = U*G*b;

h = 1;

Nx = 1001;
Ny = 1001;
Nz = 3;

M = 1.414;

beta = sqrt(M^2 - 1);

%% perturbation velocity potential function
% doublet line form, far field

syms X Y Z

r = @(x, y, z) sqrt(x.^2 - beta.^2 .*(y.^2 + z.^2));

ymax = @(x, y, z) max(-b/2, min(b/2, y + real(sqrt((x./beta).^2 - z.^2))));

ymin = @(x, y, z) max(-b/2, min(b/2, y - real(sqrt((x./beta).^2 - z.^2))));

phi_sheet = @(x, y, z, G) real(G/(2*pi).*z./abs(z) .*(atan(x./abs(z) .* (ymax(x, y, z)-y)./r(x, ymax(x, y, z)-y, z)) - atan(x./abs(z) .* (ymin(x, y, z)-y)./r(x, ymin(x, y, z)-y, z))));

phi = @(x, y, z) 0.5.*(sign(x)+1) .*(phi_sheet(x, y, z-h, G) + phi_sheet(x, y, z+h, -G));

%% grid points

xs = linspace(-150*h, 150*h, Nx)+1e-6;
ys = linspace(-150*h, 150*h, Ny)+1e-6;
zs = linspace(-2*h, 2*h, Nz);

[xp, yp, zp] = meshgrid(xs, ys, zs);

%% potential at each point

phival = phi(xp, yp, zp);

[up, vp, wp] = gradient(phival, mean(diff(xs)), mean(diff(ys)), mean(diff(zs)));

% up = up.*(1-~phival);
% vp = vp.*(1-~phival);
% wp = wp.*(1-~phival);

%% pressure perturbation at each point
Pp = -U*up - up.^2 - vp.^2 - wp.^2;
%Pp = Pp.*(sign(Pp)+1)./2;

%% check integrated lift

Lcheck = trapz(ys, trapz(xs, reshape(Pp(:, :, 2), [Nx, Ny])))/(U*G*b);

%% make plots

figure('units', 'inches', 'papersize', [4 4])
hold on
axis equal
contourf(reshape(xp(:, :, 2)./h, [Nx, Ny]), reshape(yp(:, :, 2)./h, [Nx, Ny]), log10(abs(reshape(Pp(:, :, 2), [Nx, Ny])./(L/(2*pi*h^2)))+1e-3), -3:0.5:1, 'edgecolor', 'none')
% contourf(reshape(xp(:, :, 2)./h, [Nx, Ny]), reshape(yp(:, :, 2)./h, [Nx, Ny]), reshape(Pp(:, :, 2), [Nx, Ny])./(L/(2*pi*h^2)), linspace(-1, 1, 30), 'edgecolor', 'none')
colorbar
c = colorbar;
c.Label.String = 'Log Pressure Disturbance';
c.Label.Interpreter = 'latex';
set(c,'TickLabelInterpreter','latex')
colormap(parula(8))
caxis([-3 1])
xticks([-30:10:30])
yticks([-30:10:30])
xlabel('x/h')
ylabel('y/h')
title("Ground Pressure: M = "+M)

%print("Mach"+M+".png", '-dpng', '-r600')

% figure
% hold on
% plot(xs, reshape(Pp(((Nx+1)/2),:, 2), [Nx, 1]))
% %plot(xs, L/(2*pi*h^2).*cos(atan(xs./h)).^3)




%%


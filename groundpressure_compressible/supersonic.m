clear
close all
clc

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% input parameters


Uinf = 1e3;
G = 1;
b = 1e-3;

L = Uinf*G*b;

H = 1;
NR = 100;
NTheta = 100;
Nz = 1;

Nx = NR;
Ny = NR;

M = 2;

Beta = sqrt(M^2 - 1);

%% perturbation velocity potential function
% doublet line form, far field

syms x y z f g beta h
r = sqrt(x.^2 - beta.^2 .*(y.^2 + z.^2));

maxf = (f+g + abs(f-g))/2;
minf = (abs(g-f) -f-g)/-2;

ymax = subs(maxf, [f, g], [-b/2, subs(minf, [f, g], [b/2, y + real(sqrt((x./beta).^2 - z.^2))])]);

ymin = subs(maxf, [f, g], [-b/2, subs(minf, [f, g], [b/2, y - real(sqrt((x./beta).^2 - z.^2))])]);

phi_sheet = real(1/(2*pi).*z./abs(z) .*(atan(x./abs(z) .* (ymax-y)./subs(r, y, [ymax-y])) - atan(x./abs(z) .* (ymin-y)./subs(r, y, ymin-y))));

phi = (G*subs(phi_sheet, z, z-h) + -G*subs(phi_sheet, z, z+h));

%% grid points

Rs = logspace(-0, 2, NR)*H;
Thetas = linspace(-pi, pi, NTheta);

xs = (Rs).*cosh(Thetas').*ones(1, 1, Nz);
ys = (Rs./Beta).*sinh(Thetas').*ones(1, 1, Nz);

zs = zeros(NTheta, NR, Nz)+0*H;

% zs(:,:,1) = -2*H.*ones(NTheta, NR);
% zs(:,:,2) = 0.*ones(NTheta, NR);
% zs(:,:,3) = 2*H.*ones(NTheta, NR);

% scatter3(reshape(xs, [prod(size(xs)), 1]), reshape(ys, [prod(size(ys)), 1]), reshape(zs, [prod(size(zs)), 1]))

%% potential at each point

u = diff(phi, x);
v = diff(phi, y);
w = diff(phi, z);

ufun = matlabFunction(u);
vfun = matlabFunction(v);
wfun = matlabFunction(w);
phifun = matlabFunction(phi);

up = ufun(Beta, H, xs, ys, zs);
vp = vfun(Beta, H, xs, ys, zs);
wp = wfun(Beta, H, xs, ys, zs);
phip = phifun(Beta, H, xs, ys, zs);


% up = up.*(1-~phival);
% vp = vp.*(1-~phival);
% wp = wp.*(1-~phival);

%% pressure perturbation at each point
Pp = -Uinf*up - up.^2 - vp.^2 - wp.^2;
%Pp = Pp.*(sign(Pp)+1)./2;

%% check integrated lift

%Lcheck = trapz(ys, trapz(xs, reshape(Pp(:, :, 2), [Nx, Ny])))/(U*G*b);

%% make plots

x0 = linspace(-1, 1, 3).*max(max(max(ys)));
y0 = linspace(-1, 1, 3).*max(max(max(ys)));
[x0, y0] = meshgrid(x0, y0);

figure('units', 'inches', 'papersize', [4 4])
hold on
axis equal
contourf(x0, y0, ones(size(x0))*-3+x0*1e-9, 'edgecolor', 'none')
contourf(xs(:, :, 1)./H, ys(:, :, 1)./H, log10(abs(Pp(:, :, 1)./(L/(2*pi*H^2)))+1e-3), -3:0.01:1, 'edgecolor', 'none')
colorbar
c = colorbar;
c.Label.String = 'Log Pressure Disturbance';
c.Label.Interpreter = 'latex';
set(c,'TickLabelInterpreter','latex')
%colormap(parula(8))
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
% plot(xs, reshape(Pp(((Nx+1)/2),:, 2), [Nx, 1]))
% %plot(xs, L/(2*pi*h^2).*cos(atan(xs./h)).^3)




%%


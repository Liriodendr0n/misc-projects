clear
close all
clc

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% define velocity potential functions
syms x y z f g beta h
H = 1;
L = 1;
Uinf = 1e+3;
kz = L/Uinf;
b = 1e-3;
G = 1;

%% subsonic
phi_doublet_sub = 1/(4*pi) * z./(y.^2+z.^2).*(1+x./(sqrt(x.^2+beta.^2.*(y.^2+z.^2))));
phi_sub = kz.*subs(phi_doublet_sub, z, z-h) + -kz.*subs(phi_doublet_sub, z, z+h);

u_sub = diff(phi_sub, x);
v_sub = diff(phi_sub, x);
w_sub = diff(phi_sub, x);

ufun_sub = matlabFunction(u_sub);
vfun_sub = matlabFunction(v_sub);
wfun_sub = matlabFunction(w_sub);

%% supersonic
r = sqrt(x.^2 - beta.^2 .*(y.^2 + z.^2));
maxf = (f+g + abs(f-g))/2;
minf = (abs(g-f) -f-g)/-2;
ymax = subs(maxf, [f, g], [-b/2, subs(minf, [f, g], [b/2, y + real(sqrt((x./beta).^2 - z.^2))])]);
ymin = subs(maxf, [f, g], [-b/2, subs(minf, [f, g], [b/2, y - real(sqrt((x./beta).^2 - z.^2))])]);
phi_sheet_sup = real(1/(2*pi).*z./abs(z) .*(atan(x./abs(z) .* (ymax-y)./subs(r, y, [ymax-y])) - atan(x./abs(z) .* (ymin-y)./subs(r, y, ymin-y))));
phi_sup = (G*subs(phi_sheet_sup, z, z-h) + -G*subs(phi_sheet_sup, z, z+h));

u_sup = diff(phi_sup, x);
v_sup = diff(phi_sup, y);
w_sup = diff(phi_sup, z);

ufun_sup = matlabFunction(u_sup);
vfun_sup = matlabFunction(v_sup);
wfun_sup = matlabFunction(w_sup);

%%

%% GIANT LOOP to make the movie

%Ms = 0:0.1:0.9;

Ms = [0:0.01:0.99, 1.001, 1.01:0.01:2.24]+1e-6;

figure('name', 'vidmaker', 'units', 'pixels', 'papersize', [1280 1024])

v = VideoWriter('pressdist_test.mp4', 'MPEG-4');
v.FrameRate = 60;
open(v);

for Mind = 1:length(Ms)
    
    if Ms(Mind) < 1
        %% inputs
        NR = 100;
        NTheta = 100;
        Nz = 3;

        Nx = NR;
        Ny = NR;

        M = Ms(Mind);
        Beta = sqrt(1 - M^2);
        
        %% grid points
        Rs = logspace(-3, 2, NR)*H;
        Thetas = linspace(-pi, pi, NTheta)+1e-6;

        xs = (Rs).*cos(Thetas').*ones(1, 1, 3);
        ys = (Rs./Beta).*sin(Thetas').*ones(1, 1, 3);

        zs = ones(NTheta, NR, Nz)+1e-6;

        zs(:,:,1) = -2*H.*ones(NTheta, NR);
        zs(:,:,2) = 0.*ones(NTheta, NR);
        zs(:,:,3) = 2*H.*ones(NTheta, NR);

        %% perturbation velocities and pressures
        
        us = ufun_sub(Beta, H, xs, ys, zs);
        vs = ufun_sub(Beta, H, xs, ys, zs);
        ws = ufun_sub(Beta, H, xs, ys, zs);
        
        Pp = -Uinf*us - us.^2 - vs.^2 - ws.^2;
        
        %%
    elseif Ms(Mind) > 1
        %% inputs
        NR = 100;
        NTheta = 100;
        Nz = 3;

        Nx = NR;
        Ny = NR;

        M = Ms(Mind);

        Beta = sqrt(M^2 - 1);
        
        %% grid points
        
        Rs = logspace(-1, 2, NR)*H;
        Thetas = linspace(-pi, pi, NTheta);

        xs = (Rs).*cosh(Thetas').*ones(1, 1, 3);
        ys = (Rs./Beta).*sinh(Thetas').*ones(1, 1, 3);

        zs = ones(NTheta, NR, Nz)+1e-6;

        zs(:,:,1) = -2*H.*ones(NTheta, NR);
        zs(:,:,2) = 0.*ones(NTheta, NR);
        zs(:,:,3) = 2*H.*ones(NTheta, NR);

        %% perturbation velocities and pressures
        
        us = ufun_sup(Beta, H, xs, ys, zs);
        vs = ufun_sup(Beta, H, xs, ys, zs);
        ws = ufun_sup(Beta, H, xs, ys, zs);
        
        Pp = -Uinf*us - us.^2 - vs.^2 - ws.^2;
    end
    
    %% draw the frame
    % fill background
    axis equal
    hold on
    x0 = linspace(-1, 1, 3).*max(max(max(ys)));
    y0 = linspace(-1, 1, 3).*max(max(max(ys)));
    [x0, y0] = meshgrid(x0, y0);
    contourf(x0, y0, ones(size(x0))*-3+x0*1e-9, 'edgecolor', 'none')
    
    contourf(xs(:, :, 2)./H, ys(:, :, 2)./H, log10(abs(Pp(:, :, 2)./(L/(2*pi*H^2)))+1e-3), -3:0.01:1, 'edgecolor', 'none')
    %contourf(reshape(xp(:, :, 2)./h, [Nx, Ny]), reshape(yp(:, :, 2)./h, [Nx, Ny]), reshape(Pp(:, :, 2), [Nx, Ny])./(L/(2*pi*h^2)), linspace(-1, 1, 30), 'edgecolor', 'none')
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
    grid on
    hold off
    
    ax = gca; 
    ax.FontSize = 24;
    
    f = gcf;
    f.Position = [800 20 1024 1024];
    frame = getframe(gcf);
    writeVideo(v,frame);
    disp(Ms(Mind))
    
    
    
    clf
    
end

close(v);
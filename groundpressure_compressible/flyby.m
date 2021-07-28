clear
close all
clc

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% define velocity potential functions
syms x y z f g beta h
H = 0.1;
L = 1;
Uinf = 1e+6;
kz = L/Uinf;
b = 1e-6;
G = 1;

%% subsonic
phi_doublet_sub = 1/(4*pi) * z./(y.^2+z.^2).*(1+x./(sqrt(x.^2+beta.^2.*(y.^2+z.^2))));
phi_sub = kz.*subs(phi_doublet_sub, z, z-h) + -kz.*subs(phi_doublet_sub, z, z+h);

u_sub = diff(phi_sub, x);
% v_sub = diff(phi_sub, x);
% w_sub = diff(phi_sub, x);

% ufun_sub = matlabFunction(u_sub);
% vfun_sub = matlabFunction(v_sub);
% wfun_sub = matlabFunction(w_sub);

% Pfun_sub = matlabFunction(-Uinf*u_sub - u_sub.^2 - v_sub.^2 - w_sub.^2);
Pfun_sub = matlabFunction(-Uinf*u_sub);

%% supersonic
r = sqrt(x.^2 - beta.^2 .*(y.^2 + z.^2));
maxf = (f+g + abs(f-g))/2;
minf = (abs(g-f) -f-g)/-2;
ymax = subs(maxf, [f, g], [-b/2, subs(minf, [f, g], [b/2, y + real(sqrt((x./beta).^2 - z.^2))])]);
ymin = subs(maxf, [f, g], [-b/2, subs(minf, [f, g], [b/2, y - real(sqrt((x./beta).^2 - z.^2))])]);
phi_sheet_sup = real(1/(2*pi).*z./abs(z) .*(atan(x./abs(z) .* (ymax-y)./subs(r, y, [ymax-y])) - atan(x./abs(z) .* (ymin-y)./subs(r, y, ymin-y))));
phi_sup = (1/2 + sign(x)/2)*(G*subs(phi_sheet_sup, z, z-h) + -G*subs(phi_sheet_sup, z, z+h));

u_sup = diff(phi_sup, x);
% v_sup = diff(phi_sup, y);
% w_sup = diff(phi_sup, z);

% ufun_sup = matlabFunction(u_sup);
% vfun_sup = matlabFunction(v_sup);
% wfun_sup = matlabFunction(w_sup);

% Pfun_sup = matlabFunction(-Uinf*u_sup - u_sup.^2 - v_sup.^2 - w_sup.^2);
Pfun_sup = matlabFunction(-Uinf*u_sup);

%% "Airplane" Speeds and Locations

Ntimes = 180;

Mps = [0.31, 0.87, 0.98, 1.03, 1.41, 2.24];
Yps = linspace(-1, 1, length(Mps))*50*H;
Xps = linspace(1.5, -1.5, Ntimes)'.*Mps;

%% ground grid

Nxs = 500;
Nys = 500;

xgs = linspace(-1, 1, Nxs).*80*H;
ygs = linspace(-1, 1, Nys).*80*H;
zgs = 0;
[xgs, ygs, zgs] = meshgrid(xgs, ygs, zgs);

%% Compute pressures

P = zeros(Ntimes, length(Mps), Nxs, Nys);
% this makes it like 30 times faster
P = gpuArray(P);
for t = 1:Ntimes
    for m = 1:length(Mps)
        if Mps(m) < 1
            Beta = sqrt(1 - Mps(m)^2);
            P(t, m, :, :) = Pfun_sub(Beta, H, xgs - Xps(t, m), ygs - Yps(m), zgs);
        else
            Beta = sqrt(Mps(m)^2 - 1);
            P(t, m, :, :) = Pfun_sup(Beta, H, xgs - Xps(t, m), ygs - Yps(m), zgs);
        end
        % disp(m)
    end
    disp(t)
end
% sum over all the aircraft
P = sum(P, 2);
P = gather(P);

%% test plot

% contourf(xgs, ygs, log10(abs(reshape(P(1,:,:), [Nxs, Nys])./(L/(2*pi*H^2)))+1e-3), -3:0.01:1, 'edgecolor', 'none')

%% animation

figure('name', 'vidmaker', 'units', 'pixels', 'papersize', [1280 1024]/2)

v = VideoWriter('flyby_test.mp4', 'MPEG-4');
v.FrameRate = 60;
open(v);

for ti = 1:Ntimes
    axis equal
    contourf(xgs, ygs, log10(abs(reshape(P(ti,:,:), [Nxs, Nys])./(L/(2*pi*H^2)))+1e-3), -3:0.01:1, 'edgecolor', 'none')
    colorbar
    c = colorbar;
    c.Label.Interpreter = 'latex';
    set(c,'TickLabelInterpreter','latex')
    %colormap(parula(8))
    caxis([-3 1])
    xticks([-8:2:8])
    yticks([-8:2:8])
    xlim([-8, 8])
    ylim([-8, 8])
    
    ax = gca; 
    ax.FontSize = 16;
    
    f = gcf;
    f.Position = [800 20 640 480];
    frame = getframe(gcf);
    writeVideo(v,frame);
    clf
end

close(v);

%%
%log10(abs(reshape(P(1,:,:), [Nxs, Nys])./(L/(2*pi*H^2)))+1e-3)

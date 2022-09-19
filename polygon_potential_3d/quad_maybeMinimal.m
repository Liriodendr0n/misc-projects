clear
close all
clc

%%

C = [1 1 -1; -1 1 -1; -1 -1 -1; 1 -1 -1; 1 -1 1; -1 -1 1; -1 1 1; 1 1 1]';
%C = [1 1 1; 1 -1 -1; -1 -1 1; -1 1 -1]';

% C = C + 1e-3*randn(size(C));

% grid
Nx = 10;
Ny = 10;
Nz = 10;

d = 0.001;

xg = linspace(-2, 2, Nx);
yg = linspace(-2, 2, Ny);
zg = linspace(-2, 2, Nz);

xgs = xg;
ygs = yg;
zgs = zg;

[xg, yg, zg] = ndgrid(xg, yg, zg);

xp = reshape(xg, [1,Nx*Ny*Nz]);
yp = reshape(yg, [1,Nx*Ny*Nz]);
zp = reshape(zg, [1,Nx*Ny*Nz]);

xyzp = [xp; yp; zp];

% central difference stencil
xyzppx = [xp+d/2; yp; zp];
xyzppy = [xp; yp+d/2; zp];
xyzppz = [xp; yp; zp+d/2];
xyzpmx = [xp-d/2; yp; zp];
xyzpmy = [xp; yp-d/2; zp];
xyzpmz = [xp; yp; zp-d/2];


% potential evaluation and differentiation
Csh = circshift(C, 1, 2);

uvwp = zeros(3,Nx*Ny*Nz);
uvwb = zeros(3,Nx*Ny*Nz);
for i = 1:Nx*Ny*Nz
    o(i) = pgonSolidAngle(C, xyzp(:,i));
    Opx = pgonSolidAngle(C, xyzppx(:,i));
    Opy = pgonSolidAngle(C, xyzppy(:,i));
    Opz = pgonSolidAngle(C, xyzppz(:,i));
    Omx = pgonSolidAngle(C, xyzpmx(:,i));
    Omy = pgonSolidAngle(C, xyzpmy(:,i));
    Omz = pgonSolidAngle(C, xyzpmz(:,i));
    uvwp(:,i) = [(Opx-Omx)/d; (Opy-Omy)/d; (Opz-Omz)/d]/(4*pi);
end

up = uvwp(1,:);
vp = uvwp(2,:);
wp = uvwp(3,:);

ub = uvwb(1,:);
vb = uvwb(2,:);
wb = uvwb(3,:);

upn = up./vecnorm(uvwp, 2, 1);
vpn = vp./vecnorm(uvwp, 2, 1);
wpn = wp./vecnorm(uvwp, 2, 1);

ug = reshape(up, [Nx, Ny, Nz]);
vg = reshape(vp, [Nx, Ny, Nz]);
wg = reshape(wp, [Nx, Ny, Nz]);

O = reshape(o, [Nx, Ny, Nz]);

% why permute?????
ug = permute(ug, [2,1,3]);
vg = permute(vg, [2,1,3]);
wg = permute(wg, [2,1,3]);
O = permute(O, [2,1,3]);

% streamline starts

Nxs = 5;
Nys = 5;
Nzs = 1;

[xs, ys, zs] = ndgrid(linspace(-2, 2, Nxs), linspace(-2, 2, Nys), linspace(-2, 2, Nzs));

%% plots

colors = lines(7);

figure
hold on
axis equal
grid on
grid minor

view(3)
% view([1 1 1])
xlim([-2, 2])
ylim([-2, 2])
zlim([-2, 2])
plot3([C(1,:),C(1,1)], [C(2,:),C(2,1)], [C(3,:),C(3,1)])
% quiver3(xp, yp, zp, upn, vpn, wpn)
% surf(xg, yg, O/1000, 'edgecolor', 'none')
% scatter3(xp, yp, zp, [], o, 'marker', '.')

% splt = streamline(xgs, ygs, zgs, ug, vg, wg, xs, ys, zs, [0.01, 300]);
% for i = 1:Nxs*Nys*Nzs
%     splt(i).Color = colors(2,:);
% end

% sslc1 = streamslice(xgs, ygs, zgs, ug, vg, wg, [], 0, [], 0.2);
% for i = 1:length(sslc1)
%     sslc1(i).Color = colors(2,:);
% end




% print('cubeVortStreamSlice.png', '-dpng', '-r300')

% minimal surface???? (processed for plot smoothness)
fimplicit3(@(x,y,z) mod(pgonSolidAngle(C, [x; y; z]), 4*pi)-2*pi, [-2, 2], 'edgecolor', 'none')
% xlim([-20, 20])
% ylim([-20, 20])
% zlim([-20, 20])



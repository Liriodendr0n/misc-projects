clear
close all
clc

%%

% C = [1 1 -1; -1 1 -1; -1 -1 -1; 1 -1 -1; 1 -1 1; -1 -1 1; -1 1 1; 1 1 1]';
C = [1 1 1; 1 -1 -1; -1 -1 1; -1 1 -1]';

Nx = 10;
Ny = 10;
Nz = 10;

d = 0.001;

xg = linspace(-2, 2, Nx);
yg = linspace(-2, 2, Ny);
zg = linspace(-2, 2, Nz);

[xg, yg, zg] = ndgrid(xg, yg, zg);

xp = reshape(xg, [1,Nx*Ny*Nz]);
yp = reshape(yg, [1,Nx*Ny*Nz]);
zp = reshape(zg, [1,Nx*Ny*Nz]);

xyzp = [xp; yp; zp];

xyzpdx = [xp+d; yp; zp];
xyzpdy = [xp; yp+d; zp];
xyzpdz = [xp; yp; zp+d];

uvwp = zeros(3,Nx*Ny*Nz);
for i = 1:Nx*Ny*Nz
    O = curveSolidAngle(C, xyzp(:,i));
    Odx = curveSolidAngle(C, xyzpdx(:,i));
    Ody = curveSolidAngle(C, xyzpdy(:,i));
    Odz = curveSolidAngle(C, xyzpdz(:,i));
    uvwp(:,i) = [(Odx-O)/d; (Ody-O)/d; (Odz-O)/d];
end
up = uvwp(1,:);
vp = uvwp(2,:);
wp = uvwp(3,:);

upn = up./vecnorm(uvwp, 2, 1);
vpn = vp./vecnorm(uvwp, 2, 1);
wpn = wp./vecnorm(uvwp, 2, 1);


figure
hold on
axis equal

view(3)
plot3([C(1,:),C(1,1)], [C(2,:),C(2,1)], [C(3,:),C(3,1)])
% scatter3(xp, yp, zp)
% quiver3(xp, yp, zp, upn, vpn, wpn)

fimplicit3(@(x,y,z) curveSolidAngle(C, [x; y; z])+6.2, [-1, 1], 'meshdensity', 100, 'edgecolor', 'none')
% xlim([-20, 20])
% ylim([-20, 20])
% zlim([-20, 20])



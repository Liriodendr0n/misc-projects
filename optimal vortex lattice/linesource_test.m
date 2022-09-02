clear
close all
clc

%%

G = @(x, y, z) 1./sqrt(x.^2 + y.^2 + z.^2)*1/(2);

phi = @(x, y, z) integral(@(t) G(x-t, y, z), 0, 2, 'arrayValued', true);


N = 30;

h = 1e-6;

Z0 = 0;

axlim = 3;

[xgrid, ygrid, zgrid] = ndgrid(linspace(-axlim, axlim, N), Z0, linspace(-axlim, axlim, N));

xx = reshape(xgrid, [1, N^2]);
yy = reshape(ygrid, [1, N^2]);
zz = reshape(zgrid, [1, N^2]);

uG = @(x, y, z) (phi(x-h/2, y, z) - phi(x+h/2, y, z))/h;
vG = @(x, y, z) (phi(x, y-h/2, z) - phi(x, y+h/2, z))/h;
wG = @(x, y, z) (phi(x, y, z-h/2) - phi(x, y, z+h/2))/h;

uvwG = [uG(xx, yy, zz); vG(xx, yy, zz); wG(xx, yy, zz)];

uvwL = iLineSource([xx; yy; zz], [0; 0; 0], repmat([0; 1; 0], [1, N^2]), 2*pi, 0);
% uvwBS = cross(uvwBS, repmat([2; 0; 0], [1, N^2]), 1);

uL = uvwL(1,:);
vL = uvwL(2,:);
wL = uvwL(3,:);



figure
hold on
grid on
grid minor
axis equal
view(3)
xlim([-axlim axlim])
ylim([-axlim axlim])
zlim([-axlim axlim])
% plot3([-1, 1], [0, 0], [0, 0])
% quiver3(xx, yy, zz, uG(xx, yy, zz), vG(xx, yy, zz), wG(xx, yy, zz), 0)
quiver3(xx, yy, zz, uL, vL, wL, 0)
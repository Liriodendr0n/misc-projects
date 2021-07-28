clear
close all
clc

AR = 3;
S = 1;
lambda = 0.2;
Lambda = deg2rad(33);
Gamma = deg2rad(0);

N = 30;

xyz = [0; 0; 0];

[xyzV, xyzC, xyzD, cCs, cVs] = wingGeom(xyz, N, AR, S, lambda, Lambda, Gamma);

figure
hold on
grid on
axis equal
% bound vortex
plot3(xyzV(1,:), xyzV(2,:), xyzV(3,:), 'color', [0, 0.4470, 0.7410])
for i = 1:N+1
    plot3([xyzV(1,i), 10], [xyzV(2,i), xyzV(2,i)], [xyzV(3,i), xyzV(3,i)], 'color', [0, 0.4470, 0.7410])
end
% control points
scatter3(xyzC(1,:), xyzC(2,:), xyzC(3,:), 'marker', 'x')
% downwash points
scatter3(xyzD(1,:), xyzD(2,:), xyzD(3,:),  'marker', '.')
plot3(xyzV(1,:)-cVs/4, xyzV(2,:), xyzV(3,:), 'color', 'black', 'linewidth', 2)
plot3(xyzV(1,:)+3*cVs/4, xyzV(2,:), xyzV(3,:), 'color', 'black', 'linewidth', 2)
plot3([xyzV(1,1)-cVs(1)/4, xyzV(1,1)+3*cVs(1)/4], [xyzV(2,1), xyzV(2,1)], [xyzV(3,1), xyzV(3,1)], 'color', 'black', 'linewidth', 2)
plot3([xyzV(1,end)-cVs(end)/4, xyzV(1,end)+3*cVs(end)/4], [xyzV(2,end), xyzV(2,end)], [xyzV(3,end), xyzV(3,end)], 'color', 'black', 'linewidth', 2)
xlim([-sqrt(AR*S)/2, 3*sqrt(AR*S)/2])
view(3)
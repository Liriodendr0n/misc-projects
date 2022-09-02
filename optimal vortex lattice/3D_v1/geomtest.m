clear
close all
clc

%%

Nc = 8;
Ns = 4;

[xiV, wVc, xiC, wCc] = chordLattice(Nc, true);
[etaV, wVs, etaC, wCs] = spanLattice(Ns, true, true);

@(xi, eta) [N1(xi, eta), N2(xi, eta), N3(xi, eta), N4(xi, eta)];

xyz2 = [3; 5; 0];
xyz3 = [1; -5; 0];
xyz4 = [-1; -5; 0];
xyz1 = [-3; 5; 0];

xyz = [xyz1, xyz2, xyz3, xyz4];

[xyzL, xyzR, xyzC] = buildLattice(xyz, xiV, xiC, etaV, etaC);

figure
hold on
axis equal
view(3)
grid on
grid minor
scatter3(xyzC(1,:), xyzC(2,:), xyzC(3,:))
scatter3(xyzL(1,:), xyzL(2,:), xyzL(3,:))
scatter3(xyzR(1,:), xyzR(2,:), xyzR(3,:))
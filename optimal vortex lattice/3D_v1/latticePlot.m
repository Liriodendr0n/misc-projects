clear
close all
clc

%%

Nch = 5;
Nsp = 25;

alpha = deg2rad(5);
beta = deg2rad(0);

[xiV, ~, xiC, ~] = chordLattice(Nch, true);
[etaV, wV, etaC, wC] = spanLattice(Nsp, true, true);

% panel corner points defined counter clockwise from front right (-xi +eta)
cornerPts = [-6 3 0; 0 -3 0; 1 -3 0; -5 3 0]';

[xyzL, xyzR, xyzC, xyzTL, xyzTR, ns] = buildLattice(cornerPts, xiV, xiC, etaV, etaC);

A = vorlatA(xyzC, xyzL, xyzR, ns);
b = vorlatRHS(alpha, ns);

G = A\b;


figure
hold on
grid on
grid minor
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
view(3)

% bound vortices
plot3([xyzTL(1,:); xyzL(1,:); xyzR(1,:); xyzTR(1,:)], ...
      [xyzTL(2,:); xyzL(2,:); xyzR(2,:); xyzTR(2,:)], ...
      [xyzTL(3,:); xyzL(3,:); xyzR(3,:); xyzTR(3,:)])

% collocation points
plot3(xyzC(1,:), xyzC(2,:), xyzC(3,:), 'marker', '.', 'linestyle', 'none')

% circulations
quiver3(xyzC(1,:), xyzC(2,:), xyzC(3,:), zeros(1,Nch*Nsp), zeros(1,Nch*Nsp), G')
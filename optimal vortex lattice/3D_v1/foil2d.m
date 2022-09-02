clear
close all
clc

%% Geometry

% angle of attack
alpha = asin(0.5/(2*pi));

vinf = [cos(alpha); 0; sin(alpha)];

% sweep
Lambda = deg2rad(0);

% chordwise vortices
m = 2;

% leading edge
LEV = true;

% camberline coeffs
b = [0.1, 0.1];

% % thickness coeffs
% t = 0*[0.00 -0.5 -0.*0 -0.1*0];

%% Lattice

[xiV, wV, xiC, wC] = chordLattice(m, LEV);

xiV = xiV(:)';
xiC = xiC(:)';

%% Camberline

% flattened camberline
zflag = 0;

% basis functions
nCk = @(n, k) gamma(n+1)./(gamma(k+1).*gamma(n-k+1));
B = @(x, n, k) nCk(n, k) .* (1+x).^k .* (1-x).^(n-k);
C = @(x, c1, c2) (1+x).^(c1) .* (1-x).^(c2);% ./ ((c1./(c1+c2)).^c1 .* (c2./(c1+c2)).^c2);

% camberline weights
c = [1, 1];

%% Discretization

% camberline
zeta = @(xi) b*B(xi, length(b)-1, (0:length(b)-1)').*C(xi, c(1), c(2));
dz = @(xi) (zeta(xi+0.5e-2)-zeta(xi-0.5e-2))./(1e-2);
theta = @(xi) atan(dz(xi));

dzVs = dz(xiV);
dzCs = dz(xiC);

thetaVs = theta(xiV);
thetaCs = theta(xiC);

ns = [-dzCs; zeros(size(xiC)); ones(size(xiC))];
nVs = [-dzVs; zeros(size(xiV)); ones(size(xiV))];
ts = [ones(size(xiC)); zeros(size(xiC)); dzCs];

%% Solution

xyzV = [xiV; zeros(1, m); zflag*zeta(xiV)];
xyzC = [xiC; zeros(1, m); zflag*zeta(xiC)];
bhat = [sin(Lambda); cos(Lambda); 0].*ones(3, m);

A = tatA(xyzC, xyzV, bhat, ns);

RHS = tatRHS(alpha, ns);

G = A\RHS;

%% Post Processing

if LEV == false
    % for case without leading edge lumped vortex, interpolated (expensive)
    gLE = interpMatrix(xiV, -1, 0, 0)*(G./wV)*pi;
elseif LEV == true
    % with leading edge lumped vortex, direct (cheap)
    gLE = m*G(1);
end


N = sum(G)*cos(alpha);
T = -dzVs*(G)*cos(alpha) - cos(Lambda)^-1 * (gLE)^2/(2*pi);
DL = [cos(alpha), 0, sin(alpha); 0, 0, 0; -sin(alpha), 0, cos(alpha)]*[T; 0; N];

L = DL(3);
D = DL(1);


DL
% sum(G.*(1+cV));
% sum(G)

%% Plots

figure
hold on
grid on
grid minor
axis equal
xlabel('\xi')
ylabel('\zeta')
fplot(zeta, [-1, 1])
% fplot(dz, [0, 1])
scatter(xyzV(1,:), zflag*xyzV(3,:), 'marker', 'o')
scatter(xyzC(1,:), zflag*xyzC(3,:), 'marker', 'x')

% direct K-Z forces
% quiver(xyzV(1,:), zflag*xyzV(3,:), KZ(1,:), KZ(3,:))

% linearized surface pressures
% quiver(xyzV(1,:)', zflag*xyzV(3,:)', -dzVs'.*(G.*(1+cV))*cos(alpha), (G.*(1+cV))*cos(alpha))
% quiver(xyzV(1,1), zflag*xyzV(3,1), -(m*G(1))^2/pi, 0)
quiver(xyzV(1,:)', zflag*xyzV(3,:)', -m*dzVs'.*G*cos(alpha), m*G*cos(alpha), 0)
quiver(-1, 0, -m*(m*G(1))^2/pi, 0)

% plot(xiC, wC, 'marker', '.', 'linestyle', 'none')
% plot(xiV, wV, 'marker', '.', 'linestyle', 'none')

% plot(xiV, G'./wiV.*sqrt((1-xiV)./xiV))

legend('Camberline', 'Vortices', 'Control Points', 'Surface Pressures', 'Leading Edge Suction', 'location', 'best')

% ylim([-0.5, 1.5])

% plot([0, xiC], [0, cumsum(S')/2+zeta(xiC)], 'color', 'black')
% plot([0, xiC], [0, -cumsum(S')/2+zeta(xiC)], 'color', 'black')

% figure
% hold on
% plot(xiV, Sigma)
% plot(xiV, cV)

%% Influence Matrix
function [A] = tatA(xyzC, xyzV, bhat, n)

N = size(xyzC, 2);
% t = [zeros(1, N); ones(1, N); zeros(1, N)];

[I, J] = find(ones(N));

a = iBiotSavart(xyzC(:,I), xyzV(:,J), bhat(:,J), 1, eps);

A1 = dot((a), n(:,I), 1);
A = reshape(A1, N, N);

end

%% Boundary Condition Vector

function [RHS] = tatRHS(alpha, n)

vinf = [cos(alpha), 0, sin(alpha)]';

N = size(n, 2);


RHS1 = ones(N,1);
J = find(RHS1);

RHS = dot(repmat(vinf, 1, N), n(:,J), 1)';

end

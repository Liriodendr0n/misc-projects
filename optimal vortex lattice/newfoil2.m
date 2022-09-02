clear
close all
clc

%% Geometry

% angle of attack
alpha = asin(0/(2*pi));

vinf = [cos(alpha); 0; sin(alpha)];
uinf = vinf(1);
winf = vinf(3);

% chordwise vortices
m = 10;

% camberline coeffs
b = [0.1]/(pi);

% thickness coeffs (sadly unstable when simultaneous)
t = 0*[0.1];

%% Lattice

% Quadrature points (vortices and sources)
xiQ = cos((m:-1:1)/m * pi);
wiQ = [0.5, ones(1,m-1)].*(1-xiQ)*pi/m;

% collocation points
xiC = cos((2*(m-1:-1:0)+1)/(2*m) * pi);
wiC = (1-xiC)*pi/m;

%% thickness and camber

nCk = @(n, k) gamma(n+1)./(gamma(k+1).*gamma(n-k+1));
B = @(x, n, k) nCk(n, k)./(2.^n) .* (1+x).^k .* (1-x).^(n-k);

C = @(x, c1, c2) (1+x).^(c1) .* (1-x).^(c2);

t = t/1.5^1.5/0.5^0.5*2;

tfun = @(x) t*B(x, length(t)-1, (0:length(t)-1)').*C(x, 0.5, 1.5);
zfun = @(x) b*B(x, length(b)-1, (0:length(b)-1)').*C(x, 1, 1);

h = 1e-4;
dtfun = @(x) 1/h * (-3/2*tfun(x) + 2*tfun(x+h) - 1/2*tfun(x+2*h));
dzfun = @(x) 1/h * (zfun(x+h/2) - zfun(x-h/2));

Ufun = @(x) zfun(x) + tfun(x)/2;
Lfun = @(x) zfun(x) - tfun(x)/2;

dUfun = @(x) 1/h * (-3/2*Ufun(x) + 2*Ufun(x+h) - 1/2*Ufun(x+2*h));
dLfun = @(x) 1/h * (-3/2*Lfun(x) + 2*Lfun(x+h) - 1/2*Lfun(x+2*h));

nU = [-dUfun(xiC); ones(1,m)];
nL = [dLfun(xiC); -ones(1,m)];

nG = [-dzfun(xiC); zeros(1,m); ones(1,m)];
nS = [-dtfun(xiC)./sqrt((1-xiC)./(1+xiC)); zeros(1,m); zeros(1,m)];

V = 1./(xiC'-xiQ)/(2*pi);

bG = sin(alpha)+cos(alpha)*nG(1,:)';
% bS = nS(1,:)';

% G = V\b;
% G = V\bG;
G = V\bG;
% G = real(X);
% S = imag(X);

N = sum(G)*cos(alpha);
T = -dzfun(xiQ)*(G)*cos(alpha) - (m*G(1))^2/(2*pi);
DL = [cos(alpha), 0, sin(alpha); 0, 0, 0; -sin(alpha), 0, cos(alpha)]*[T; 0; N]

figure
hold on
axis equal
grid on
grid minor

plot(xiQ, G)
% plot(xiQ, S)

% fplot(Ufun, [-1, 1])
% fplot(zfun, [-1, 1])
% fplot(Lfun, [-1, 1])
% quiver(xiC, Ufun(xiC), nU(1,:), nU(2,:))
% quiver(xiC, Lfun(xiC), nL(1,:), nL(2,:))
% quiver(xiC, zeros(1,m), nU(1,:)-nL(1,:), nU(2,:)-nL(2,:))
% quiver(xiC, zeros(1,m), nU(1,:)+nL(1,:), nU(2,:)+nL(2,:))
% quiver(xiC, zeros(1,m), -dzfun(xiC), ones(1,m))
% quiver(xiC, zeros(1,m), -dtfun(xiC), zeros(1,m))



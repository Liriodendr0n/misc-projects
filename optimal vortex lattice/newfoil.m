clear
close all
clc

%% Geometry

% angle of attack
alpha = asin(1/(2*pi));

vinf = [cos(alpha); 0; sin(alpha)];
uinf = vinf(1);
winf = vinf(3);

% chordwise vortices
m = 4;

% camberline coeffs
b = [0];

% thickness coeffs (sadly unstable when simultaneous)
t = [0.4];

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
twpoly = @(x) t*B(x, length(t)-1, (0:length(t)-1)');
zfun = @(x) b*B(x, length(b)-1, (0:length(b)-1)').*C(x, 1, 1);

h = 1e-6;
dtfun = @(x) 1/h * (-3/2*tfun(x) + 2*tfun(x+h) - 1/2*tfun(x+2*h));
dzfun = @(x) 1/h * (zfun(x+h/2) - zfun(x-h/2));

W = @(x) sqrt((1-x)./(1+x));
Z = @(x) 1-x.^2;
T = @(x) (1-x).^1.5.*(1+x).^0.5;

zq = zfun(xiQ);
tq = tfun(xiQ);
dzq = dzfun(xiQ);
dtq = dtfun(xiQ);

zc = zfun(xiC);
tc = tfun(xiC);
dzc = dzfun(xiC);
dtc = dtfun(xiC);

%% influence matrices

Wc = diag(W(xiC));

Tc = tfun(xiC)';
Zc = zfun(xiC)';

Tcq = interpMatrix(xiQ, xiC, 0.5, -0.5);
TqQ = diag(1./wiQ);

[~, Dc_t] = interpMatrix(xiC, xiC+eps, 1.5, 0.5);
[~, Dc_z] = interpMatrix(xiC, xiC+eps, 1, 1);

[~, D_ut] = interpMatrix(xiC, xiC+eps, 1.5, 0.5);
[~, D_gz] = interpMatrix(xiC, xiC+eps, 1.5, 0.5);
[~, D_gt] = interpMatrix(xiC, xiC+eps, -1, 0);
[~, D_uz] = interpMatrix(xiC, xiC+eps, 1, 1);



AcQ = 1./(xiC'-xiQ)/(2*pi);

A_w = AcQ;
A_l = Tcq*TqQ;

A_ut = D_ut*(Tc./T(xiC)'.*AcQ);
A_gz = D_gz*Zc./T(xiC)'.*Tcq*TqQ;

A_gt = D_gt*(Tc.*Tcq*TqQ);

A_uz = D_uz*(Zc.*AcQ);


% A = [A_l - 0*A_ut, -0*A_gz; -0*A_uz, A_w + 1/4*A_gt];
% RHS = [uinf*Dc_t*(tc./T(xiC))'; winf - uinf*Dc_z*(zc./Z(xiC))'];
% 
% X = A\RHS;
% 
% L = X(1:m);
% G = X(m+1:2*m);


A = A_w + 1/4*A_gt;
RHS = winf - uinf*Dc_z*(zc./Z(xiC))';

G = A\RHS;

% L;
% sum(L);
sum(G)

%%

% dzVs = dz(xiQ);
% dzCs = dz(xiC);

figure
hold on
axis equal

fplot(@(x) tfun(x)/2 + zfun(x), [-1, 1], 'color', 'black')
fplot(@(x) -tfun(x)/2 + zfun(x), [-1, 1], 'color', 'black')
fplot(@(x) zfun(x), [-1, 1], 'color', 'black')

% scatter(xiC, A_gt*G)
scatter(xiQ, m*G)
% scatter(xiQ, m*L)

% scatter(xiC, A_gt*G)

% scatter(xiC, 10*Dc_t*(Tc.*AcQ)*L)
% scatter(xiC, Tc)
% scatter(xiC, Zc)
% scatter(xiC, Dc_t*Zc./T(xiC').*Tcq*TqQ*G)






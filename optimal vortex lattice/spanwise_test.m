clear
close all
clc

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% planform shape


alpha = deg2rad(1);
% alpha = asin(1/(2*pi));
beta = 0;

vinf = [cos(alpha); -sin(beta); sin(alpha)];

AR = 5;
S = 1;
lambda = 1;
Lambda = deg2rad(0);
Gamma = deg2rad(0);

% % convair delta
% AR = 2.23;
% S = 1;
% lambda = 0;
% Lambda = deg2rad(53.5);
% Gamma = deg2rad(0);

% % warren - 12
% AR = 2*sqrt(2);
% S = 2*sqrt(2);
% lambda = 1/3;
% Lambda = deg2rad(49.64);
% Gamma = deg2rad(0);

% chordwise
m = 2;

% spanwise
n = 8;

%% shape functions

b = sqrt(AR*S);
cMAC = S/b;
cfun = @(eta) 2*S/((1+lambda)*b) * (1-(1-lambda)*abs(eta));
xLEfun = @(eta) 1/4*cfun(0) - 1/4*cfun(eta) + tan(Lambda)*abs(b/2 * eta);
xTEfun = @(eta) 1/4*cfun(0) + 3/4*cfun(eta) + tan(Lambda)*abs(b/2 * eta);

xWfun = @(eta, xi) (1-xi).*xLEfun(eta) + xi.*xTEfun(eta);
yWfun = @(eta, xi) b/2 * eta;
zWfun = @(eta, xi) tan(Gamma)*abs(b/2 * eta);

% Ts = [cos(alpha), 0, sin(alpha); 0, 1, 0; -sin(alpha), 0, cos(alpha)];

%% spanwise lattice

% % ideal for smooth wings (rectangle and delta)
% etaV = cos((2*(1:n+1)-1)/(2*n+2) * pi);
% etaC = cos((1:n)/(n+1) * pi);


% ideal for kinked wings (swept or strongly tapered)
[etaV, etaC] = spanTip(floor(n/2));
etaV = [flip(etaV)', -etaV(2:end)'];
etaC = [flip(etaC)', -etaC'];

%% chordwise lattice

% % vortices
% xiV = sin((2*(m:-1:1) - 1)/(2*m+1) *pi/2).^2;
% % controls
% xiC = sin((2*(m:-1:1))/(2*m+1) *pi/2).^2;

% vortices (radau, accurate LE suction)
xiV = (1+cos((1:m)/m * pi))/2;
% controls
xiC = (1+cos((2*(0:m-1)+1)/(2*m) * pi))/2;

%% vortex lattice

% control points
[etaC, xiC] = meshgrid(etaC, xiC);
% left vortex endpoints
[etaVL, xiVL] = meshgrid(etaV(1:end-1), xiV);
% right vortex endpoints
[etaVR, xiVR] = meshgrid(etaV(2:end), xiV);
% vortex "mid"points
[etaVM, xiVM] = meshgrid(etaC(1,:), xiV);

xC = reshape(xWfun(etaC, xiC)', 1, m*n);
yC = reshape(yWfun(etaC, xiC)', 1, m*n);
zC = reshape(zWfun(etaC, xiC)', 1, m*n);

xVL = reshape(xWfun(etaVL, xiVL)', 1, m*n);
yVL = reshape(yWfun(etaVL, xiVL)', 1, m*n);
zVL = reshape(zWfun(etaVL, xiVL)', 1, m*n);

xVR = reshape(xWfun(etaVR, xiVR)', 1, m*n);
yVR = reshape(yWfun(etaVR, xiVR)', 1, m*n);
zVR = reshape(zWfun(etaVR, xiVR)', 1, m*n);

xVM = reshape(xWfun(etaVM, xiVM)', 1, m*n);
yVM = reshape(yWfun(etaVM, xiVM)', 1, m*n);
zVM = reshape(zWfun(etaVM, xiVM)', 1, m*n);

xyzC = [xC; yC; zC];
xyzVL = [xVL; yVL; zVL];
xyzVR = [xVR; yVR; zVR];
xyzVM = [xVM; yVM; zVM];

Lgamma = xyzVL - xyzVR;

ns = cross((xyzVL-xyzC), (xyzVR-xyzC))./vecnorm(cross((xyzVL-xyzC), (xyzVR-xyzC)));

ts = cross(xyzVR-xyzVL, ns, 1)./vecnorm(cross(xyzVR-xyzVL, ns, 1), 2, 1);
cs = cfun(etaC);

xhat = [ones(1, n*m); zeros(1, n*m); zeros(1, n*m)];
yhat = [zeros(1, n*m); ones(1, n*m); zeros(1, n*m)];
zhat = [zeros(1, n*m); zeros(1, n*m); ones(1, n*m)];


% aileron
% TEind = n*(m-1):n*m;
% 
% ns(1, TEind(10:floor(n/2)-5)) = sind(20);
% ns(3, TEind(10:floor(n/2)-5)) = cosd(20);

%% linear system

% local normal velocity
A = vorlatA(xyzC, xyzVL, xyzVR, ns);

% trefftz plane normalwash
Atfz = vorlatAtfz(xyzC, xyzVL, xyzVR, ns);

% local chordwise perturbation velocity
Acdw = vorlatAcdw(xyzC, xyzVL, xyzVR, ns);

% neumann B.C.
RHS = vorlatRHS(alpha, ns, zeros(1, m*n)');

% forward velocity evaluation

Vx = vorlatA(xyzVM, xyzVL, xyzVR, xhat);
Vy = vorlatA(xyzVM, xyzVL, xyzVR, yhat);
Vz = vorlatA(xyzVM, xyzVL, xyzVR, zhat);

G = A\RHS;
w = Atfz*G;

Gnm = reshape(G, [n, m]);

G1 = Gnm(:,m);

Gn = sum(Gnm, 2);

uvw = vinf - [Vx*G, Vy*G, Vz*G]';

KZ = cross(uvw, xyzVL-xyzVR, 1).*G';
kz = cross(uvw, (xyzVL-xyzVR)./vecnorm((xyzVL-xyzVR), 2, 1), 1).*G';


% trefftz plane

% wWeights = pi/(n+1) * sin((1:n)/(n+1)*pi).^2;

%% force calculation

% trefftz forces (far field)

CL1 = 2*sum((yVL-yVR).*G')/S
CD1 = sum((yVL-yVR).*w'.*G')/S;
e1 = CL1^2/(pi*AR*CD1);

% local forces (near field)

N = S*CL1/2 * cos(alpha);
LEsucc = (( m*G1 ).^2./dot(ts(:,end-n+1:end),[-1; 0; 0].*ones(3, n))'./(cs(1,:)'*pi));
T = (yVR(1:n)-yVL(1:n))*LEsucc;

DL = [cos(alpha), sin(alpha); -sin(alpha), cos(alpha)]*[T; N];

CL2 = 2*DL(2)/S;
CD2 = 2*DL(1)/S;
e2 = CL2^2/(pi*AR*CD2);

DLKZ = [cos(alpha), 0, sin(alpha); 0, 0, 0; -sin(alpha), 0, cos(alpha)]*sum(KZ, 2);

CL3 = 2*DLKZ(3)/S;
CD3 = 2*DLKZ(1)/S;
e3 = CL3^2/(pi*AR*CD3);

%% trefftz plot

colors = lines(7);

figure%('name', 'spanplot', 'units', 'inches', 'papersize', [8, 6], 'paperposition', [0 0 8 6])
hold on
grid on
grid minor

plot(yC(1:n), sum(reshape(G, [n, m]), 2))
plot(yC(1:n), sum(reshape(-w/m, [n, m]), 2))
plot(yC(1:n), 10*LEsucc)
plot(yC(1:n), -10*sum(reshape(w.*G, [n, m]), 2))
plot(yC(1:n), -20*sum(reshape(G, [n, m]), 2)*sin(alpha) + 20*LEsucc)
legend('$C_l$', '$w$', '$s\cdot10$', '$C_{d, tfz}\cdot10$', '$C_{d, les}\cdot10$', 'location', 'best')
% legend('$C_l$', '$C_{d, \mathrm{tfz}}\cdot10$', '$C_{d, \mathrm{les}}\cdot10$', 'location', 'northeast')

% title('Warren 12 Planform')

% print('spanplot.pdf', '-dpdf', '-painters')

%% plot geometry

% figure%('name', 'vlm', 'units', 'inches', 'papersize', [8, 6])
% hold on
% grid on
% grid minor
% axis equal
% set(gca, 'ydir', 'reverse')
% 
% plot3([b*ones(1, m*n); xVL; xVR; b*ones(1, m*n)], [yVL; yVL; yVR; yVR], [zVL; zVL; zVR; zVR], 'color', 'black')
% plot3(xC, yC, zC, 'marker', '.', 'color', 'black', 'linestyle', 'none')
% 
% % % surface pressure
% quiver3(xVM, yVM, zVM, ns(1,:).*G', ns(2,:).*G', ns(3,:).*G', 2, 'color', colors(1,:))
% % plot3(reshape(xyzVM(1,:), [n, m]), reshape(xyzVM(2,:), [n, m]), m*[ones(1, m-1), 2].*Gnm./(1-xiVM')./sqrt(1-etaVM'.^2), 'color', colors(1,:))
% % plot3(reshape(xyzVM(1,:), [n, m])', reshape(xyzVM(2,:), [n, m])', m*[ones(1, m-1), 2]'.*Gnm'./(1-xiVM)./sqrt(1-etaVM.^2), 'color', colors(1,:))
% 
% % % LE suction
% quiver3(xVM(end-n+1:end), yVM(end-n+1:end), zVM(end-n+1:end), ...
%      ts(1,end-n+1:end).*LEsucc', ts(2,end-n+1:end).*LEsucc', ts(3,end-n+1:end).*LEsucc', 0.5, 'color', colors(2,:))
% % plot3(reshape(xyzVM(1,end-n+1:end)+m.*G(end-n+1:end)'.*ts(1,end-n+1:end), [n/2, 2]), reshape(xyzVM(2,end-n+1:end)+m.*G(end-n+1:end)'.*ts(2,end-n+1:end), [n/2, 2]), zeros(2, n/2), 'color', colors(2,:))
% % plot3((xyzVM(1,end-n+1:end)+[zeros(1, n); m.*G(end-n+1:end)'.*ts(1,end-n+1:end)]), (xyzVM(2,end-n+1:end)+[zeros(1, n); m.*G(end-n+1:end)'.*ts(2,end-n+1:end)]), zeros(2, n)', 'color', colors(2,:))
% 
% 
% 
% % K-Z Forces
% % quiver3(xVM, yVM, zVM, kz(1,:), kz(2,:), kz(3,:), 2, 'color', colors(1,:))
% 
% 
% plot3(xWfun([0, 1, 1, 0, -1, -1, 0], [0, 0, 1, 1, 1, 0, 0]), ...
%       yWfun([0, 1, 1, 0, -1, -1, 0], [0, 0, 1, 1, 1, 0, 0]), ...
%       zWfun([0, 1, 1, 0, -1, -1, 0], [0, 0, 1, 1, 1, 0, 0]), 'color', 'black', 'linestyle', '--')
% 
% view(3)
% xlabel('x')
% ylabel('y')
% zlabel('z')
% % set(gca, 'ydir', 'reverse')
% % view([0, 0, 1])
% % print('vlm.pdf', '-dpdf', '-painters')

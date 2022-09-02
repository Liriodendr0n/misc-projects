clear
close all
clc

%%

% alpha = asin(0.1/(2*pi));
% alpha = 2/pi;
% alpha = 0.1;
alpha = deg2rad(1);


n = 30;
c = 1;
b = 10;

AR = b/c;

% ideal for smooth wings (rectangle and delta)
etaQ = cos((2*(1:n+1)-1)/(2*n+2) * pi);
etaC = cos((1:n)/(n+1) * pi);
etaM = etaQ(1:n)+diff(etaQ)/2;

wQ = pi/(n+1) .*sqrt(1-etaQ.^2);
wC = pi/(n+1) *sin((1:n)/(n+1)*pi).^2 ./sqrt(1-etaC.^2);

% cs = sqrt(1-etaC.^2)';
cs = 1*ones(n, 1);
alpha0s = 0*ones(n, 1);

% [M, DM, IM] = interpMatrix(etaQ, etaC, -0.5, -0.5);

% classical lifting line
W = diff(-1./(etaC'-etaQ)/(4*pi), 1, 2)*(2/b);
A = eye(n) - pi*cs.*W;
RHS = pi*cs.*(alpha0s+alpha);
G = A\RHS;
L = wC*G
D = wC*(G.*(-W*G));
% L = diff(-etaQ)*G
% D = diff(etaQ)*(G.*(W*G));
e0 = L^2/(pi*AR*D);

% % modified
% W = 1./(etaC'-etaQ)/(4*pi)*2/b;
% A = [tril(ones(n,n+1)) - pi*cs.*W; ones(1,n+1)];
% % A = [-IM*(n+1)/pi - pi*cs.*W; ones(1,n+1)];
% RHS = [pi*cs.*(alpha0s+alpha); 0];
% g = A\RHS;
% G = cumsum(g);
% G = G(1:n);
% L = wC*G
% D = wC*(G.*(-W*g));
% % L = diff(-etaQ)*G
% % D = diff(etaQ)*(G.*(W*g));
% e0 = L^2/(pi*AR*D);

figure
hold on
grid on
grid minor
plot(etaC, 2*G)
% plot(etaC, W*g)
plot(etaC, W*G)
% plot(etaC, -2*(n+1)/pi*IM*g)

clear
close all
clc

%%

% exponentially convergent or even exact for nice geometries

alpha = asin(1/(2*pi));

N = 100;

[etaQ, wQ, etaC, wC] = spanLattice(N, true, true);

% cfun = @(eta) sqrt(1 - eta.^2);
cfun = @(eta) 1 - 0*abs(eta);

b = 10;

cs = cfun(etaC);


% shed vorticity influence matrix
A = -1./(etaQ'-etaC) * 1/(2*pi) * 2/b * pi/(N+1);

% bound to shed interpolation
[~, B] = interpMatrix(etaC, etaQ, 0.5, 0.5);
B = B.*sqrt(1-etaQ.^2)./sqrt(1-etaC'.^2);

G = (A*B + eye(N)./(pi*cs))\(ones(N,1)*alpha);

CL = wC'*(G./sqrt(1-etaC.^2))

figure
hold on
plot(etaC, G)
plot(etaC, -A*B*G)
% plot(etaC, G./sqrt(1-etaC.^2))

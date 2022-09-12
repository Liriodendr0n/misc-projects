clear
close all
clc

%%

% exponentially convergent or even exact for nice geometries

alpha = asin(1/(2*pi));

N = 10;

Ltip = false;
Rtip = false;

[etaQ, wQ, etaC, wC] = spanLattice(N, Ltip, Rtip);

WQ = (1-etaQ).^(0.5*Rtip).*(1+etaQ).^(0.5*Ltip);
WC = (1-etaC).^(0.5*Rtip).*(1+etaC).^(0.5*Ltip);

% cfun = @(eta) sqrt(1 - eta.^2);
cfun = @(eta) 1 - 0.0*abs(eta);

% k = 0.1;
% cfun = @(eta) 1 - 0.5*(1+sqrt(eta.^2+k^2)-sqrt(1-k^2));

afun = @(eta) 0.00*erf(1e3*(eta-0.5));

b = 10;

cs = cfun(etaC);
as = afun(etaC);

% cs([ceil(N/2)]) = cs([ceil(N/2)])+0.00;

% shed vorticity influence matrix
A = -1./(etaQ'-etaC) * 1/(2*pi) * 2/b .* wQ';

% so much faster (don't use 2, ill-conditioned)
[~, B] = interpMatrix3(etaC, etaQ, 0.5*Rtip, 0.5*Ltip);

[muTip] = interpMatrix3(etaC, [-1, 1], 0.5*Rtip, 0.5*Ltip);
B(1,:) = B(1,:) + muTip(1,:)./wQ(1);
B(end,:) = B(end,:) - muTip(end,:)./wQ(end);

B = B.*WQ./WC';

% LLT linear system solution
G = (A*B + eye(N)./(pi*cs))\(as + alpha);

% G = G./sqrt(1-etaC'.^2)';

% CL quadrature rule
CL = (wC./WC)'*G
CD = (wC./WC)'*(G.*A*B*G)


figure
hold on
% plot(etaC, G)
plot(etaC, G./cs)
plot(etaC, -A*B*G)
% plot(etaC, G./sqrt(1-etaC.^2))


% cInt = @(x) interpMatrix3(etaC, x, 0.5, 0.5)*(cs./sqrt(1-etaC.^2));
% figure
% hold on
% fplot(cfun, [-1, 1])
% fplot(cInt, [-1, 1])

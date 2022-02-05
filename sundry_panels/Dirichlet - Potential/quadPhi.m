function [Phi] = quadPhi(x0, y0, x1, y1, x2, y2, beta1, beta2, mu0, mu1, mu2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% cubic hermite panels (makes paneling so much easier)
% up to quadratic singularity strength

%% translate and rotate

panel_ang = -atan2(y2-y1, x2-x1);
panel_len = sqrt((x2-x1).^2 + (y2-y1).^2);

x12 = 0.5*x1 + 0.5*x2;
y12 = 0.5*y1 + 0.5*y2;

x = 2*(cos(panel_ang).*x0 - sin(panel_ang).*y0 - x12.*cos(panel_ang)+y12.*sin(panel_ang))/panel_len;
y = 2*(sin(panel_ang).*x0 + cos(panel_ang).*y0 - x12.*sin(panel_ang)-y12.*cos(panel_ang))/panel_len;

%% kernels

phiDx = @(xi, eta) -xi./(xi.^2 + eta.^2);
phiDy = @(xi, eta) -eta./(xi.^2 + eta.^2);

zeta = @(t) beta1/4 * (t+1).*(1-t).^2 + beta2/4 * (t+1).^2.*(t-1);

dzeta = @(t) (t-1).^2/4+(t+1).^2/4+((t*2+2).*(t-1))/4+(t*2-2).*(t/4+1/4);

ds = @(t) sqrt(1 + dzeta(t).^2);

mu = @(t) mu0 + t.*mu1 + t.^2.*mu2;

%% integral

integrand = @(xi, eta, t) mu(t).*ds(t) .* (cos(atan(dzeta(t))).*phiDy(xi-t, eta-zeta(t)) - sin(atan(dzeta(t))).*phiDx(xi-t, eta-zeta(t)));

Phi = 1/(2*pi) * integral(@(t) integrand(x, y, t), -1, 1);


end


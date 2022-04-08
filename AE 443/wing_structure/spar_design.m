clear
close all
clc

%% properties

% airfoil section thickness
tcr = 0.17;
tcs = 0.17;
tct = 0.13;

% wing geometry
S = 500;
AR = 4.0;
lambda = 0.75;
kappa = 0.1;

b = sqrt(AR*S);
c0 = 2*b/(AR*(kappa + lambda - kappa*lambda + 1));


cfun = @(y) interp1([-b/2, -kappa*b/2, kappa*b/2, b/2], [lambda*c0, c0, c0, lambda*c0], y);
tcfun = @(y) interp1([-b/2, -kappa*b/2, 0, kappa*b/2, b/2], [tct, tcs, tcr, tcs, tct], y);

ys = linspace(-b/2, b/2, 31);

%% aero loads

wing.AR = AR;
wing.S = S;
wing.lambda = lambda;
wing.kappa = kappa;

flight.v_inf = 800;
flight.rho_inf = 0.00238;
flight.alpha = deg2rad(5-3);
flight.omega = deg2rad(0);

section.a0Lflap = deg2rad(2.2);
section.a0LailR = deg2rad(2.2);
section.a0LailL = deg2rad(2.2);

[Lprime] = LLT(wing, section, flight);
centerind = ceil(31/2);
Lprimecalc = Lprime(1:centerind);
yscalc = ys(1:centerind);

Ltot = trapz(ys, Lprime);
Wplane = 14000;
n_test = Ltot/Wplane;

% spar loads
sparShear = -cumtrapz(sign(yscalc).*yscalc, Lprimecalc);
sparBend = cumtrapz(sign(yscalc).*yscalc, sparShear);

% rib loads
ribload = @(x) 2/pi * (1-x)./(sqrt(x.*(1-x)));
Nribs = 15;
ysrib = linspace(-b/2, 0, Nribs);
Dy = b./(2*Nribs);

%LExs = linspace(-cfun(ysrib)*0.25, 0, 10);
%TExs = linspace(cfun(ysrib)*0.75, 0, 10);

%ribLoadLE = Lprime*ribload(LExs + cfun(ysrib))

%ribShear = 
%ribBend = 

%% spar geometry

sigma = 5760000/1.15 % fatigue strength, 40 ksi into psf
rho_spar = 168.55549;   % lbm/ft^3

% main spar shear web height, a function of thickness and chord
hweb = cfun(yscalc)'.*tcfun(yscalc)';
% shear web thickness
tweb = sparShear./(sigma*hweb);
Wweb = trapz(yscalc, hweb.*tweb*rho_spar);
Aweb = hweb.*tweb;

% spar cap area
Fcap = abs(sparBend)./hweb;
Acap = Fcap/(sigma);
Wcaps = 2.*trapz(yscalc, Acap.*rho_spar);

lcap = hweb./2;
tcap = Acap./lcap;
% spar weight
Wspar = 2*(Wweb + Wcaps)

%% rib geometry

% rib height constant and equal to shear web height for simplicity


%hrib = cfun(ysrib)'.*tcfun(ysrib)';
%trweb = 

%Wspar

Ltot/4600

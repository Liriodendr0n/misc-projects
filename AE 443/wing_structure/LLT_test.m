clear
close all
clc

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% problem data

v_inf = 1;        % ft/s
rho_inf = 1;
omega_roll = 0.0;     % rad/s

S = 1;        % ft^2
AR = 6.0;

% tip taper ratio
lambda = 1;
%shoulder position
kappa = 0.001;

% wing angle of attack
alpha_geo = deg2rad(9.2);
% section zero lift angles of attack (section slope 2pi/rad)
alpha0Lflap = deg2rad(0);
alpha0LailR = deg2rad(0);
alpha0LailL = deg2rad(0);

b = sqrt(AR*S);
c0 = 2*b/(AR*(kappa + lambda - kappa*lambda + 1));
cfun = @(y) c0.*(heaviside(y + kappa*b/2).*heaviside(kappa*b/2 - y) + ...
    heaviside(y - kappa*b/2).*(y*(2*(lambda-1))/(b*(1-kappa)) - (kappa*(lambda-1))/(1-kappa) + 1) + ...
    heaviside(-y - kappa*b/2).*(y*(2*(1-lambda))/(b*(1-kappa)) + (kappa*(1-lambda))/(1-kappa) + 1));

alpha0Lfun = @(y) alpha0Lflap.*(heaviside(y + kappa*b/2).*heaviside(kappa*b/2 - y)) + ...
    alpha0LailR.*heaviside(y - kappa*b/2) + alpha0LailL.*heaviside(-y - kappa*b/2);

CLopt = 0.6;
Oopt = (2*(1+lambda)*CLopt)/(pi*(2*pi));
twistfun = @(y) -(Oopt.* (1 - sqrt(1-(2*y./b).^2)./(cfun(y)./c0)));

N = 100;

ellip = true;

%% set up lifting line geometry

Dy = b/N;                               % span of each horseshoe vortex
ys = linspace(-b/2+Dy/2, b/2-Dy/2, N);  % y position of each control point
cs = cfun(ys);
alpha_aero = alpha_geo - alpha0Lfun(ys) - 0*twistfun(ys)- omega_roll*ys/v_inf; % alpha_geo - alpha_ind - alpha_0L

%% set up LHS vector b

for i = 1:N
    LHS(i,:) = pi*v_inf*cs(i)*alpha_aero(i);
end

%% set up influence matrix A

for i = 1:N
    for j = 1:N
        
        % downwash from each horseshoe vortex: left and right shed vortices
        % effect of horseshoe j on horseshoe i
        k(i,j) = 1/(4*pi*(ys(j)-ys(i) + Dy/2)) - 1/(4*pi*(ys(j)-ys(i) - Dy/2));
        
        % build the A matrix
        if i == j
            A(i,j) = (1 + pi*cs(i)*k(i,j));
        else
            A(i, j) = pi*cs(i)*k(i,j);
        end
        
    end
end

%% solve linear system A*Gs = LHS
% vector of circulation strengths of horseshoe vortices
tic
Gs = A\LHS;
toc

%% calculate downwash
alpha_ind = -(Gs./(v_inf*cs'*pi) - alpha_aero');
w = v_inf*alpha_ind;

%% calculate aero coefficients
Cl = 2*Gs./(v_inf*cs');
CL = sum(2*Gs*Dy)./(v_inf*S);
Lprime = 0.5*rho_inf*v_inf^2*cs'.*Cl;

Cdi = Cl.*alpha_ind;
CDi = sum(cs'.*Cdi.*Dy)/S;

Cd0 = 0.005;
Cd2 = 0.00005;
Cdv = Cd0 + Cd2*rad2deg((alpha_aero'-alpha_ind)).^2;
Cd = Cdi + Cdv;
CDv = sum(cs'.*Cdv.*Dy)/S;
CD = CDi + CDv;

Cm0 = 0;
Cm1 = 0;
Cm = Cm0 + Cm1*rad2deg(alpha_aero'-alpha_ind);
CM = sum(cs'.*Cm.*Dy)/S;

%% calculate aero ratios
LDi = CL/CDi;
LD = CL/CD;

%% print results to command window
disp("CL = " + CL)
disp("CDi = " + CDi)
disp("CDv = " + CDv)
disp("L/Di = " + LDi)
disp("L/D = " + LD)

%% plot results
figure
hold on
% plot(ys, Gs)
plot(ys, Cl)
ylim([0 1])



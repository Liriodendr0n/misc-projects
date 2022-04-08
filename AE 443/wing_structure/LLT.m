function [Lprime, cs, ys, alpha0L] = LLT(wing, section, flight)
%% use

% FLIGHT --> alpha rho_inf v_inf omega
% WING --> AR S lambda kappa
% SECTION --> a0Lflap a0LailR a0LailL Cd0 Cd2 Cm0 Cm1

%% problem data

v_inf = flight.v_inf;        % ft/s
omega_roll = flight.omega;     % rad/s
rho_inf = flight.rho_inf;       % slug/ft^3

S = wing.S;        % ft^2
AR = wing.AR;

% tip taper ratio
lambda = wing.lambda;
%shoulder position
kappa = wing.kappa;

% wing angle of attack
alpha_geo = flight.alpha;
% section zero lift angles of attack (section slope 2pi/rad)
alpha0Lflap = section.a0Lflap;
alpha0LailR = section.a0LailR;
alpha0LailL = section.a0LailL;

b = sqrt(AR*S);
c0 = 2*b/(AR*(kappa + lambda - kappa*lambda + 1));
cfun = @(y) c0.*(heaviside(y + kappa*b/2).*heaviside(kappa*b/2 - y) + ...
    heaviside(y - kappa*b/2).*(y*(2*(lambda-1))/(b*(1-kappa)) - (kappa*(lambda-1))/(1-kappa) + 1) + ...
    heaviside(-y - kappa*b/2).*(y*(2*(1-lambda))/(b*(1-kappa)) + (kappa*(1-lambda))/(1-kappa) + 1));

alpha0Lfun = @(y) alpha0Lflap.*(heaviside(y + kappa*b/2).*heaviside(kappa*b/2 - y)) + ...
    alpha0LailR.*heaviside(y - kappa*b/2) + alpha0LailL.*heaviside(-y - kappa*b/2);

CLopt = 0.6;
lambda = 0.45;
Oopt = (2*(1+lambda)*CLopt)/(pi*(2*pi));
twistfun = @(y) -(Oopt.* (1 - sqrt(1-(2*y./b).^2)./(cfun(y)./c0)));

N = 31;

%% set up lifting line geometry

Dy = b/N;                               % span of each horseshoe vortex
ys = linspace(-b/2+Dy/2, b/2-Dy/2, N)';  % y position of each control point
cs = cfun(ys);
alpha_aero = alpha_geo - alpha0Lfun(ys') - twistfun(ys') - omega_roll*ys'/v_inf; % alpha_geo - alpha_ind - alpha_0L
alpha0L = alpha0Lfun(ys');

%% set up LHS vector

for i = 1:N
    LHS(i,:) = pi*v_inf*cs(i)'*alpha_aero(i);
end

%% set up influence matrix A

for i = 1:N
    for j = 1:N
        
        % downwash from each horseshoe vortex: left and right shed vortices
        % effect of horseshoe j on horseshoe i
        k(i,j) = 1/(4*pi*(ys(j)'-ys(i)' + Dy/2)) - 1/(4*pi*(ys(j)'-ys(i)' - Dy/2));
        
        % build the A matrix
        if i == j
            A(i,j) = (1 + pi*cs(i)'*k(i,j));
        else
            A(i, j) = pi*cs(i)'*k(i,j);
        end
        
    end
end

%% solve linear system A*Gs = LHS
% vector of circulation strengths of horseshoe vortices
Gs = A\LHS;

%% calculate downwash
alpha_ind = -(Gs./(v_inf*cs*pi) - alpha_aero');
w = v_inf*alpha_ind;

%% calculate aero coefficients
Cl = 2*Gs./(v_inf*cs);
CL = sum(2*Gs*Dy)./(v_inf*S);
Lprime = rho_inf*v_inf*Gs;

% Cdi = Cl.*alpha_ind;
% CDi = sum(cs'.*Cdi.*Dy)/S;
% 
% Cd0 = section.Cd0;
% Cd2 = section.Cd2;
% Cdv = Cd0 + Cd2*rad2deg((alpha_aero'-alpha_ind)).^2;
% Cd = Cdi + Cdv;
% CDv = sum(cs'.*Cdv.*Dy)/S;
% CD = CDi + CDv;
% 
% Cm0 = section.Cm0;
% Cm1 = section.Cm1;
% Cm = Cm0 + Cm1*rad2deg(alpha_aero'-alpha_ind);
% CM = sum(cs'.*Cm.*Dy)/S;


end


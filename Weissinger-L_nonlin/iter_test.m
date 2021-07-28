clear
close all
clc

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

load naca2412_1e5.mat

Clfun = @(alpha) interp1(deg2rad(naca24121e5.alpha), naca24121e5.CL, alpha);
Cdfun = @(alpha) interp1(deg2rad(naca24121e5.alpha), naca24121e5.CD, alpha);
% Clfun = @(alpha) 2*pi*sin(alpha) + 0.48;
% Cdfun = @(alpha) sin(alpha).^2 + 0.006;

%% geometry and conditions
tic
% conditions

alpha = deg2rad(5);

% geometry
AR = 6.5;
S = 1;
lambda = 0.5;
Lambda = deg2rad(0);
Gamma = deg2rad(0);

N = 100;

xyz = [0; 0; 0];

b = sqrt(AR*S);
c = S/b;

% error tolerance (iterate until Delta CL smaller than this)
CLtol = 1e-15;

% build wing geometry
[xyzV, xyzC, xyzD, cCs, cVs, ns] = wingGeom(xyz, N, AR, S, lambda, Lambda, Gamma);

%% compute influence matrices
% control point influence matrix
A = weisslA(xyzC, xyzV, ns);

% trefftz plane downwash matrix
Atfz = weisslAtrefftz(xyzD, xyzV, ns);

t_build = toc;

%% initial linear solution
tic
% initial RHS
RHS = weisslRHS(alpha, ns, 0*ones(N,1));
% solve for circulation strengths
G = A\RHS;
% compute initial Cl and spanload
Cl = 2*G./cCs;
Ccl = 2*G/c;
% compute induced angle of attack
alpha_ind = -atan2(0 - Cl/(2*pi), cos(alpha)) - alpha;
% compute transpiration velocity
Vtransp = sin(alpha + alpha_ind) - Clfun(alpha + alpha_ind)/(2*pi);

t_init = toc;

%% iterative nonlinear solution

i = 2;

err = 1;

tic
while err > CLtol
    % main loop
    % build RHS from flight state and previous transpiration
    RHS(:,i) = weisslRHS(alpha, ns, Vtransp(:,i-1));
    % solve for circulation strengths
    G(:,i) = A\RHS(:,i);
    % predict Cls and spanload from circulation strengths
    Cl(:,i) = 2*G(:,i)./(cCs);
    Ccl(:,i) = 2*G(:,i)/c;
    % compute induced angle from current predicted Cl and previous Vt
    alpha_ind(:,i) = -atan2(-Vtransp(:,i-1) - Cl(:,i)/(2*pi), cos(alpha)) - alpha;
    % correct Cls from current induced angle
    Cl(:,i) = 0.5*Clfun(alpha + alpha_ind(:,i)) + 0.5*Clfun(alpha + alpha_ind(:,i-1));
    Ccl(:,i) = Cl(:,i).*cCs/c;
    % compute current Vt
    Vtransp(:,i) = sin(alpha + alpha_ind(:,i)) - Clfun(alpha + alpha_ind(:,i))/(2*pi);
    
    % integrate Cls along span
    CL(1,i) = sum(diff(xyzV(2,:)'/b).* Ccl(:,i));
    % relative CL error (current - previous)
    err(:,i) = abs(CL(1,i) - CL(1,i-1));
    
    % update index
    i = i+1;
end

% trefftz downwash
w = -0.5*Atfz*G(:,end);
% induced drag and dragload
Cdi = w.*Cl(:,end);
Ccdi(:,i) = w.*Ccl(:,end);

% integrate full wing drag
CDi = trapz(xyzC(2,:)/b, Cdi(:,end).*cCs./c);
CDv = trapz(xyzC(2,:)/b, Cdfun(alpha + alpha_ind(:,end)).*cCs./c);
CD = CDi + CDv;
e = (CL(:,end).^2)/(pi*AR*CDi)


t_soln = toc;
iterend = i-2;

%% results

disp("build time = " + round(t_build*1000, 2) + "ms")
% disp("init time = " + round(t_init*1000, 2) + "ms")
disp("solution time = " + round(t_soln*1000, 2) + "ms")
disp("solution iterations = " + iterend)
disp("CL = " + CL(end))
disp("CDi = " + CDi)
disp("CD = " + CD)
disp("L/Di = " + CL(end)/CDi)
disp("L/D = " + CL(end)/CD)

%% figure summary

figure('name', 'results', 'units', 'inches', 'papersize', [10 10], 'paperposition', [0 0 10 10 ])

sgtitle("$\alpha = \,$" + rad2deg(alpha) + "$^\circ$, " + ...
    "$AR = \,$" + AR + ", " + ...
    "$\lambda = \,$" + lambda + ", " + ...
    "$\Lambda = \,$" + rad2deg(Lambda) + "$^\circ$, " + ...
    "$\Gamma = \,$" + rad2deg(Gamma) + "$^\circ$, " + ...
    "$N = \,$" + N)

% convergence
subplot(2, 2, 1)
hold on
plot(0:i-2, err, 'marker', '.')
set(gca, 'Yscale', 'log')
grid on
xlabel('Iteration')
ylabel('$\Delta C_L$ Error')

% local Cl
subplot(2, 2, 2)
hold on
grid on
grid minor
plot(2*xyzC(2,:)/b, Cl(:,end))
plot(2*xyzC(2,:)/b, Clfun(alpha + alpha_ind(:,end)), 'linestyle', '--')
ylim([-1.0, 2.0])
xlabel('$2y/b$')
ylabel('$c_l$')
plot(2*xyzC(2,:)/b, Cl(:,end).*cCs/c, ':', 'color', 'black')
legend('Prediction', 'Correction', 'location', 'best')

% data fitting
subplot(2, 2, 3)
hold on
grid on
grid minor
fplot(@(a) Clfun(deg2rad(a)), [-90, 90])
scatter(rad2deg(alpha+alpha_ind(:,end)), Cl(:,end), 'marker', 'x')
fplot(@(a) Cdfun(deg2rad(a)), [-90, 90])
scatter(rad2deg(alpha + alpha_ind(:,end)), Cdfun(alpha + alpha_ind(:,end)), 'marker', 'x')
ylim([-1.0, 2.0])
xlim([-30, 30])
xlabel('$\alpha$')
ylabel('$c_l \, , c_d$')
legend('Section Data', 'Span Stations', 'location', 'best')

% geometry
subplot(2, 2, 4)
hold on
grid on
axis equal
% control points
scatter3(xyzC(1,:), xyzC(2,:), xyzC(3,:), 10*ones(N,1), ones(N,1).*[0.8500, 0.3250, 0.0980], 'marker', 'x')
% downwash points
%scatter3(xyzD(1,:), xyzD(2,:), xyzD(3,:), 15*ones(N,1), ones(N,1).*[0.9290, 0.6940, 0.1250],  'marker', '.')
% bound vortex
plot3(xyzV(1,:), xyzV(2,:), xyzV(3,:), 'color', [0, 0.4470, 0.7410], 'linewidth', 1.5)
for i = 1:N+1
    plot3([xyzV(1,i), 1e3], [xyzV(2,i), xyzV(2,i)], [xyzV(3,i), xyzV(3,i)], 'color', [0.3010, 0.7450, 0.9330])
end
% trailing vortices
plot3(xyzV(1,:)-cVs/4, xyzV(2,:), xyzV(3,:), 'color', 'black', 'linewidth', 1.5)
plot3(xyzV(1,:)+3*cVs/4, xyzV(2,:), xyzV(3,:), 'color', 'black', 'linewidth', 1.5)
plot3([xyzV(1,1)-cVs(1)/4, xyzV(1,1)+3*cVs(1)/4], [xyzV(2,1), xyzV(2,1)], [xyzV(3,1), xyzV(3,1)], 'color', 'black', 'linewidth', 1.5)
plot3([xyzV(1,end)-cVs(end)/4, xyzV(1,end)+3*cVs(end)/4], [xyzV(2,end), xyzV(2,end)], [xyzV(3,end), xyzV(3,end)], 'color', 'black', 'linewidth', 1.5)
xlim([-0.5*b, 1.0*b])
ylim([-0.75*b, 0.75*b])
zlim([-0.333*b, 0.667*b])
view(3)
xlabel('x')
ylabel('y')
zlabel('z')
legend('Control Points', 'Bound Vortex', 'Trailing Vortices', 'location', 'best')

%% print summary

%print('resultsVis.pdf', '-dpdf', '-painters')

%% extra figures

% figure
% hold on
% plot(2*xyzC(2,:)/b, Cl(:,end))
% plot(2*xyzC(2,:)/b, Clfun(alpha + alpha_ind(:,end)), 'linestyle', '--')



%%

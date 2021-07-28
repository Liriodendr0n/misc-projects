clear
close all
clc

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

load MSc175e6.mat

% Clfun = @(alpha) heaviside(deg2rad(25) - abs(alpha)).*interp1(deg2rad(naca24121e5.alpha), naca24121e5.CL, alpha, 'pchip')...
%     + heaviside(abs(alpha) - deg2rad(25)).*(2.1.*sin(alpha)+0.31).*cos(alpha);

Clfun = @(alpha) interp1(deg2rad(MSc175e6.alpha), MSc175e6.CL, alpha, 'pchip');

%% problem parameters
alpha = deg2rad(9.2);
% ||v_infty|| is assumed to be 1

xyz = [0; 0; 0];
N = 1000;

AR = 6.6;
S = 1;
lambda = 0.5;
Lambda = deg2rad(0);
Gamma = deg2rad(0);


b = sqrt(AR*S);
c = S/b;
%% build geometry
[xyzV, xyzC, xyzD, cCs, cVs, ns] = wingGeom(xyz, N, AR, S, lambda, Lambda, Gamma);

%% build and solve linear system

CLdes = 1;

Vtransp = 0*(sin(alpha) - CLdes/(2*pi));

RHS = weisslRHS(alpha, ns, 1*Vtransp*ones(N,1));

A = weisslA(xyzC, xyzV, ns);
Atrefftz = weisslAtrefftz(xyzD, xyzV, ns);

G = A\RHS;

% local Cl
Cl = 2*G./cCs;
% local Loading
Ccl = 2*G/c;
% perpendicular Cl
Clp = Cl/(cos(Lambda)^2);

%% downwash
% local induced angle (based on local Cl and Vt) (THIS IS THE ONE)
%alphai = -atan2(-Vtransp - Cl/(2*pi), cos(alpha)) - alpha;

% trefftz downwash
ai = atan(-0.5*Atrefftz*G);

Cdi = sin(ai).*Cl/cos(Gamma);
Ccdi = sin(ai).*Ccl/cos(Gamma);


%Vtransp1 = (sin(alpha) - Cl/(2*pi));

%% results
CL = sum(diff(xyzV(2,:))/b * Ccl)
CDi = sum(diff(xyzV(2,:))/b * Ccdi)

e = (CL^2)/(pi*AR*CDi)


%% aero figures 
figure
hold on
grid on

yyaxis left
plot(xyzC(2,:), Cl)
plot(xyzC(2,:), Ccl)
ylim([0, 1.0])

yyaxis right
plot(xyzC(2,:), Cdi)
plot(xyzC(2,:), Ccdi)
ylim([0, 1.0]/10)

legend('Section $C_l$', 'Spanload', 'location', 'best')

figure
hold on
plot(xyzC(2,:), -atan(2*ai))
%plot(xyzC(2,:), alphai)
ylim([-0.2, 0.1])

%% visualization
% figure
% hold on
% grid on
% axis equal
% % bound vortex
% plot3(xyzV(1,:), xyzV(2,:), xyzV(3,:), 'color', [0, 0.4470, 0.7410])
% for i = 1:N+1
%     plot3([xyzV(1,i), 10], [xyzV(2,i), xyzV(2,i)], [xyzV(3,i), xyzV(3,i)], 'color', [0, 0.4470, 0.7410])
% end
% % control points
% scatter3(xyzC(1,:), xyzC(2,:), xyzC(3,:), 'marker', 'x')
% % downwash points
% scatter3(xyzD(1,:), xyzD(2,:), xyzD(3,:),  'marker', '.')
% plot3(xyzV(1,:)-cVs/4, xyzV(2,:), xyzV(3,:), 'color', 'black', 'linewidth', 2)
% plot3(xyzV(1,:)+3*cVs/4, xyzV(2,:), xyzV(3,:), 'color', 'black', 'linewidth', 2)
% plot3([xyzV(1,1)-cVs(1)/4, xyzV(1,1)+3*cVs(1)/4], [xyzV(2,1), xyzV(2,1)], [xyzV(3,1), xyzV(3,1)], 'color', 'black', 'linewidth', 2)
% plot3([xyzV(1,end)-cVs(end)/4, xyzV(1,end)+3*cVs(end)/4], [xyzV(2,end), xyzV(2,end)], [xyzV(3,end), xyzV(3,end)], 'color', 'black', 'linewidth', 2)
% xlim([-sqrt(AR*S), sqrt(AR*S)])
% view(3)


%%

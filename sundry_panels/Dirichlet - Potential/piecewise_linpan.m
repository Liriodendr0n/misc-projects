clear
close all
clc

%%

% piecewise linear doublet panels, two independent unknowns per panel
% gauss collocation points


alpha = deg2rad(10);
V = 1;

Qinf = [V*cos(alpha); V*sin(alpha)];

N = 30;
coords = genNACA4([2 4 12], 1, N+1);

xp = coords(:,1);
yp = coords(:,2);

thetap = atan2(diff(-yp),-diff(xp));
ds = sqrt(diff(xp).^2 + diff(yp).^2);

xc1 = xp(1:N) + diff(xp)*(3-sqrt(3))/6; % - 1e-9*cos(thetap - pi/2);
xc2 = xp(1:N) + diff(xp)*(3+sqrt(3))/6; % - 1e-9*cos(thetap - pi/2);
yc1 = yp(1:N) + diff(yp)*(3-sqrt(3))/6; % - 1e-9*sin(thetap - pi/2);
yc2 = yp(1:N) + diff(yp)*(3+sqrt(3))/6; % - 1e-9*sin(thetap - pi/2);

xc = [xc1; xc2];
yc = [yc1; yc2];

s = [0; cumsum(sqrt(diff(xp).^2 + diff(yp).^2))];

% figure
% hold on
% axis equal
% grid on
% grid minor
% plot(xp, yp, 'marker', '.')
% plot(xc1, yc1, 'linestyle', 'none', 'marker', '.')
% plot(xc2, yc2, 'linestyle', 'none', 'marker', '.')

%% influences

% first block: constant strengths
% second block: linear strengths
% points i, singularities j

A = zeros(2*N+1);

%      -----------------------------        |mu0(1)|
%     | left points  | left points  |       | ...  |
%     | constant     | linear       |       |mu0(N)|
% A = |--------------+---------------  mu = |mu1(1)|
%     | right points | right points |       | ...  |
%     | constant     | linear       |       |mu1(N)|
%      -----------------------------        |mu0(W)|

% AICs
for i = 1:N
    for j = 1:N
        if i == j
            A(i, j) = 0.5;
            A(i+N, j) = 0.5;
            A(i, j+N) = -sqrt(3)/6;
            A(i+N, j+N) = sqrt(3)/6;
        else
            % constant strengths
            A(i, j) = polyPhi(xc1(i), yc1(i), xp(j), yp(j), xp(j+1), yp(j+1), 1, 0, 0);
            A(i+N, j) = polyPhi(xc2(i), yc2(i), xp(j), yp(j), xp(j+1), yp(j+1), 1, 0, 0);
            % linear strengths
            A(i, j+N) = polyPhi(xc1(i), yc1(i), xp(j), yp(j), xp(j+1), yp(j+1), 0, 1, 0);
            A(i+N, j+N) = polyPhi(xc2(i), yc2(i), xp(j), yp(j), xp(j+1), yp(j+1), 0, 1, 0);
        end
    end
end


% wake
for i = 1:N
    A(i, 2*N+1) = polyPhi(xc1(i), yc1(i), 1, 0, 1000, 0, 1, 0, 0);
    A(i+N, 2*N+1) = polyPhi(xc2(i), yc2(i), 1, 0, 1000, 0, 1, 0, 0);
end

% wake strength
A(2*N+1, 1) = 1; A(2*N+1, N) = -1; A(2*N+1, N+1) = 0; A(2*N+1, 2*N) = 0;
A(2*N+1, 2*N+1) = 1;

% RHS
RHS = -[[xc, yc]*Qinf; 0];

% kutta velocity
% RHS(2*N) = 0;
% A(2*N, :) = zeros(1,2*N+1);
% A(2*N, N+1) = ds(end); A(2*N, 2*N) = ds(1);


%% solve
mu = A\RHS;

% plot(xp(1:N) + diff(xp), mu(1:N))
% 
Qt = 2*mu(N+1:2*N)./ds;
Cp = 1 - Qt.^2;

figure
hold on
grid on
grid minor
plot([xp(1:N), xp(2:N+1)]', [mu(1:N)-mu(N+1:2*N), mu(1:N)+mu(N+1:2*N)]', 'linewidth', 1.5)

figure
hold on
grid on
grid minor
set(gca, 'ydir', 'reverse')
plot(xp(1:N)+0.5*diff(xp), Cp)
colors = lines(7);
% plot([xp(1:N), xp(2:N+1)]', [Cp, Cp]', 'color', colors(1,:))

Cl = -2*mu(end);
disp("Cl = " + round(Cl, 3))


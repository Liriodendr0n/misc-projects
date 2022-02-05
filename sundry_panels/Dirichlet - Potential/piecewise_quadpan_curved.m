clear
close all
clc

%%

% piecewise linear doublet panels, two independent unknowns per panel
% gauss collocation points


alpha = deg2rad(0);
V = 1;

Qinf = [V*cos(alpha); V*sin(alpha)];

N = 31;
coords = genNACA4([2 4 12], 1, N+1);

xp = coords(:,1);
yp = coords(:,2);

thetap = atan2(diff(yp), diff(xp));
dtheta = atan2(sin(diff(thetap)), cos(diff(thetap)));

ds = sqrt(diff(xp).^2 + diff(yp).^2);

np = [thetap(1); ...
      atan2((ds(1:end-1)./(ds(1:end-1)+ds(2:end))).*sin(thetap(1:end-1)) + ...
            (ds(2:end)./(ds(1:end-1)+ds(2:end))).*sin(thetap(2:end)), ...
            (ds(1:end-1)./(ds(1:end-1)+ds(2:end))).*cos(thetap(1:end-1)) + ...
            (ds(2:end)./(ds(1:end-1)+ds(2:end))).*cos(thetap(2:end))); ...
      thetap(end)] + pi/2;

beta1 = -tan((atan2(sin(thetap - np(1:end-1)+pi/2), cos(thetap - np(1:end-1)+pi/2))));
beta2 = -tan((atan2(sin(thetap - np(2:end)+pi/2), cos(thetap - np(2:end)+pi/2))));

zf = @(b1, b2, t) b1 * t.*(1-t).^2 + b2 * t.^2.*(t-1);

t1 = (1-sqrt(3/5))/2;
t2 = 1/2;
t3 = (1+sqrt(3/5))/2;

for i = 1:N
    xc1(i, 1) = xp(i) + t1*(xp(i+1)-xp(i)) + cos(thetap(i) + pi/2).*ds(i).*zf(beta1(i), beta2(i), t1)*0.99;
    xc2(i, 1) = xp(i) + t2*(xp(i+1)-xp(i)) + cos(thetap(i) + pi/2).*ds(i).*zf(beta1(i), beta2(i), t2)*0.99;
    xc3(i, 1) = xp(i) + t3*(xp(i+1)-xp(i)) + cos(thetap(i) + pi/2).*ds(i).*zf(beta1(i), beta2(i), t3)*0.99;

    yc1(i, 1) = yp(i) + t1*(yp(i+1)-yp(i)) + sin(thetap(i) + pi/2).*ds(i).*zf(beta1(i), beta2(i), t1)*0.99;
    yc2(i, 1) = yp(i) + t2*(yp(i+1)-yp(i)) + sin(thetap(i) + pi/2).*ds(i).*zf(beta1(i), beta2(i), t2)*0.99;
    yc3(i, 1) = yp(i) + t3*(yp(i+1)-yp(i)) + sin(thetap(i) + pi/2).*ds(i).*zf(beta1(i), beta2(i), t3)*0.99;
end

xc = [xc1; xc2; xc3];
yc = [yc1; yc2; yc3];

% figure
% hold on
% axis equal
% grid on
% grid minor
% plot(xp, yp, 'marker', '.')
% plot(xc1, yc1, 'linestyle', 'none', 'marker', '.')
% plot(xc2, yc2, 'linestyle', 'none', 'marker', '.')
% plot(xc3, yc3, 'linestyle', 'none', 'marker', '.')

%% influences

% first block: constant strengths
% second block: linear strengths
% points i, singularities j

A = zeros(3*N+1);

% AICs
for i = 1:N
    for j = 1:N
        % constant strengths
        A(i, j) = quadPhi(xc1(i), yc1(i), xp(j), yp(j), xp(j+1), yp(j+1), beta1(j), beta2(j), 1, 0, 0);
        A(i+N, j) = quadPhi(xc2(i), yc2(i), xp(j), yp(j), xp(j+1), yp(j+1), beta1(j), beta2(j), 1, 0, 0);
        A(i+2*N, j) = quadPhi(xc3(i), yc3(i), xp(j), yp(j), xp(j+1), yp(j+1), beta1(j), beta2(j), 1, 0, 0);
        % linear strengths
        A(i, j+N) = quadPhi(xc1(i), yc1(i), xp(j), yp(j), xp(j+1), yp(j+1), beta1(j), beta2(j), 0, 1, 0);
        A(i+N, j+N) = quadPhi(xc2(i), yc2(i), xp(j), yp(j), xp(j+1), yp(j+1), beta1(j), beta2(j), 0, 1, 0);
        A(i+2*N, j+N) = quadPhi(xc3(i), yc3(i), xp(j), yp(j), xp(j+1), yp(j+1), beta1(j), beta2(j), 0, 1, 0);
        % quadratic strengths
        A(i, j+2*N) = quadPhi(xc1(i), yc1(i), xp(j), yp(j), xp(j+1), yp(j+1), beta1(j), beta2(j), 0, 0, 1);
        A(i+N, j+2*N) = quadPhi(xc2(i), yc2(i), xp(j), yp(j), xp(j+1), yp(j+1), beta1(j), beta2(j), 0, 0, 1);
        A(i+2*N, j+2*N) = quadPhi(xc3(i), yc3(i), xp(j), yp(j), xp(j+1), yp(j+1), beta1(j), beta2(j), 0, 0, 1);
    end
end


% wake
for i = 1:N
    A(i, 3*N+1) = polyPhi(xc1(i), yc1(i), 1, 0, 1000, 0, 1, 0, 0);
    A(i+N, 3*N+1) = polyPhi(xc2(i), yc2(i), 1, 0, 1000, 0, 1, 0, 0);
    A(i+2*N, 3*N+1) = polyPhi(xc3(i), yc3(i), 1, 0, 1000, 0, 1, 0, 0);
end

% wake strength
A(3*N+1, 1) = 1; A(3*N+1, N) = -1;
% A(3*N+1, N+1) = ds(end); A(3*N+1, 2*N) = -ds(1);
% A(3*N+1, 2*N+1) = ds(end).^2; A(3*N+1, 3*N) = -ds(1).^2;

A(3*N+1, 3*N+1) = 1;

% RHS
RHS = -[[xc, yc]*Qinf; 0];

% kutta velocity
RHS(3*N) = 0;
A(3*N, :) = zeros(1,3*N+1);
A(3*N, N+1) = ds(end); A(3*N, 2*N) = ds(1);
% A(3*N, 2*N+1) = 2*ds(end); A(3*N, 3*N) = 2*ds(1);

%% solve
mu = A\RHS;

% plot(xp(1:N) + diff(xp), mu(1:N))
% 
Cl = -2*mu(end);
disp("Cl = " + round(Cl, 3))

%% plot

figure
hold on
axis equal
grid on
grid minor
for i = 1:N
    fplot(@(t) interp1(1:N+1, xp, t) + zf(beta1(i), beta2(i), t-i).*cos(thetap(i)+pi/2).*ds(i), ...
          @(t) interp1(1:N+1, yp, t) + zf(beta1(i), beta2(i), t-i).*sin(thetap(i)+pi/2).*ds(i), [i, i+1], 'linewidth', 1.5)
end
% quiver(xp, yp, -cos(np), -sin(np), 'color', 'black')
% plot(xp, yp, 'marker', '.', 'color', 'black')
plot(xc, yc, 'linestyle', 'none', 'color', 'black', 'marker', '.')

figure
plot(xp(1:N)+diff(xp)/2, mu(1:N))

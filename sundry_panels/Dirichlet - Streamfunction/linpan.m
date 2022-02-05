clear
close all
clc

%%

alpha = deg2rad(0);

Qinf = [cos(alpha); sin(alpha)];

N = 1000;
wU = [0.16 0.16 0.22 0.12];
wL = [-0.09 -0.06 -0.12 0.12];
% set to zero
dzTE = 0.00;

coords = genCST(wU, wL, dzTE, N);
% coords = genNACA4([2 4 12], 1, N);
xp = coords(:,1);
yp = coords(:,2);

ds = sqrt(diff(xp).^2 + diff(yp).^2);
% ds = sqrt(gradient(xp).^2 + gradient(yp).^2);

s = [0; cumsum(sqrt(diff(xp).^2 + diff(yp).^2))];

xc = xp;
yc = yp;

% left and right panel endpoints
xpl = xp(1:N-1);
xpr = xp(2:N);
ypl = yp(1:N-1);
ypr = yp(2:N);

% trailing edge angle
thetaTEU = -atan(3*sqrt(3)/4*wU(end)-dzTE/2);
thetaTEL = -atan(3*sqrt(3)/4*wL(end)+dzTE/2);
thetaTE = 0.5*thetaTEU + 0.5*thetaTEL;

%% AIC

A = zeros(N+1);
Psil = zeros(N, N-1);
Psir = zeros(N, N-1);
for i = 1:N
    for j = 1:N-1
        Psil(i, j) = lvs(xc(i), yc(i), xpr(j), ypr(j), xpl(j), ypl(j))*ds(j);
        Psir(i, j) = lvs(xc(i), yc(i), xpl(j), ypl(j), xpr(j), ypr(j))*ds(j);
    end
end
A(1:N, 1:N-1) = Psil;
A(1:N, 2:N) = A(1:N, 2:N) + Psir;

% trailing edge gap
if dzTE ~= 0
    for i = 1:N
        A(i, 1) = A(i, 1) - 0.5*abs(cos(thetaTE))*css(xc(i), yc(i), xp(1), yp(1), xp(N), yp(N))*dzTE;
        A(i, 1) = A(i, 1) - 0.5*abs(sin(thetaTE))*cvs(xc(i), yc(i), xp(1), yp(1), xp(N), yp(N))*dzTE;

        A(i, N) = A(i, N) + 0.5*abs(cos(thetaTE))*css(xc(i), yc(i), xp(1), yp(1), xp(N), yp(N))*dzTE;
        A(i, N) = A(i, N) + 0.5*abs(sin(thetaTE))*cvs(xc(i), yc(i), xp(1), yp(1), xp(N), yp(N))*dzTE;
    end
end

% surface streamfunction
A(:, N+1) = -1;

% kutta condition
A(N+1, :) = [1, zeros(1, N-2), 1, 0];

%% b and solution

b = [[yc, -xc]*Qinf; 0];

% TE panel extrapolation (sharp TE)

if dzTE == 0
    b(N) = 0;
    A(N, :) = [1, -2, 1, zeros(1, N-6), -1, 2, -1, 0];
end

x = A\b;

g = x(1:N);
Psi0 = x(N+1);

qt = g;
Cp = 1 - qt.^2;

Cl = 2*sum(qt.*[0.5*ds(1); 0.5*ds(1:N-2)+0.5*ds(2:N-1); ds(N-1)]);
disp("Cl = " + Cl)

%% plots

colors = lines(7);

figure
hold on
grid on
grid minor
set(gca, 'ydir', 'reverse')
plot(xp, Cp)

% figure
% hold on
% grid on
% grid minor
% plot(s/s(end), qt)

% figure
% hold on
% grid on
% grid minor
% axis equal
% fimplicit(@(x, y) linPan_PsiPlot((x-0.5)*cos(-alpha)+y*sin(-alpha)+0.5, (0.5-x)*sin(-alpha)+y*cos(-alpha), ...
%     xp, yp, g, Qinf, thetaTE, dzTE)-Psi0, [-0.5, 1.5, -1.0, 1.0], 'color', [0 0 0])
% fimplicit(@(x, y) linPan_PsiPlot((x-0.5)*cos(-alpha)+y*sin(-alpha)+0.5, (0.5-x)*sin(-alpha)+y*cos(-alpha), ...
%     xp, yp, g, Qinf, thetaTE, dzTE)-Psi0+0.01, [-0.5, 1.5, -1.0, 1.0], 'color', [0 0 0])
% plot(polyshape((xp-0.5)*cos(alpha)+yp*sin(alpha)+0.5, (0.5-xp)*sin(alpha)+yp*cos(alpha)), ...
%     'facealpha', 1, 'facecolor', [0.67 0.67 0.67])
% % xlim([-0.5, 1.5])
% xlim([0.95, 1.05])

%%



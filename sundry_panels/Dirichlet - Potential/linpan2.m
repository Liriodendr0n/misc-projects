clear
close all
clc

%%

alpha = deg2rad(0);
V = 1;

Qinf = [V*cos(alpha); V*sin(alpha)];

N = 100;
coords = genNACA4([2 4 12], 1, N);

xp = coords(:,1);
yp = coords(:,2);

s = [0; cumsum(sqrt(diff(xp).^2 + diff(yp).^2))];

thetap = atan2(diff(-yp),-diff(xp));

thetac = [atan2(0.5*sin(thetap(1))+0.5*sin(thetap(end)), 0.5*cos(thetap(1))+0.5*cos(thetap(end))); ...
          atan2(0.5*sin(thetap(1:end-1))+0.5*sin(thetap(2:end)), 0.5*cos(thetap(1:end-1))+0.5*cos(thetap(2:end))); ...
          atan2(0.5*sin(thetap(1))+0.5*sin(thetap(end)), 0.5*cos(thetap(1))+0.5*cos(thetap(end)))];

xc = xp - 1e-9*cos(thetac - pi/2);
yc = yp - 1e-9*sin(thetac - pi/2);

figure
hold on
axis equal
plot(xp, yp)
plot(xc, yc, 'linestyle', 'none', 'marker', '.')
quiver(xp, yp, cos(thetac-pi/2), sin(thetac-pi/2))

%%

Phia = zeros(N);
Phib = zeros(N);

% N panels, N+1 vertices (trailing edge doublet strength is not continuous)
% collocation points i, panels j
for i = 1:N
    for j = 1:N-1
        [Phia(i, j), Phib(i, j)] = linPhi(xc(i), yc(i), xp(j), yp(j), xp(j+1), yp(j+1), 1, 1);
    end
end

%% assemble A
A = Phia + circshift(Phib, 1, 2);


% wake
for i = 1:N
    A(i, N+1) = constPhi(xc(i), yc(i), 1, 0, 1000, 0, 1);
end

% kutta condition (TE potential jump applies to the wake)
A(N+1, :) = [1, zeros(1, N-2), -1, 1];

RHS = -[[xc, yc]*Qinf; 0];

% kutta condition (TE doublet gradient constant)
RHS(N) = 0;
A(N, :) = [-1, 1, zeros(1, N-4), -1, 1, 0];

mu = A\RHS;

% take second order central differences to smooth oscillations
Qt = diff(mu(1:N))./diff(s);
Cp = 1 - Qt.^2/(Qinf'*Qinf);


figure
hold on
% ylim([-1, 1])
% xlim([0, 1])
plot(xp(1:N-1) + 0.5*diff(xp), Cp)
set(gca, 'ydir', 'reverse')

% figure
% hold on
% plot(s, Qt)

Cl = -2*mu(end);
disp("Cl = " + round(Cl, 3))
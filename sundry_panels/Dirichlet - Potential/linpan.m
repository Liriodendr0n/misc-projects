clear
close all
clc

%%

alpha = deg2rad(10);
V = 1;

Qinf = [V*cos(alpha); V*sin(alpha)];

N = 60;
coords = genNACA4([2 4 12], 1, N);

xp = coords(:,1);
yp = coords(:,2);

s = [0; cumsum(sqrt(diff(xp).^2 + diff(yp).^2))];

thetap = atan2(diff(-yp),-diff(xp));

xc = xp(1:N-1) + 0.5*diff(xp);
yc = yp(1:N-1) + 0.5*diff(yp);

xc = xc;
yc = yc;

xc(end+1) = 1-(1e-9);
yc(end+1) = 0;

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

for i = 1:N-1
    Phia(i, i) = 0.25;
    Phib(i, i) = 0.25;
end

%% assemble A
A = Phia + circshift(Phib, 1, 2);

% for i = 1:N
%     A(i, i) = -0.25;
% end

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
Qt = gradient(mu(1:N))./gradient(s);
Cp = 1 - Qt.^2/(Qinf'*Qinf);


figure
hold on
% ylim([-1, 1])
% xlim([0, 1])
plot(xp, Cp)
grid on
grid minor
set(gca, 'ydir', 'reverse')

% figure
% hold on
% plot(s, Qt)

% mu(end)

Cl = -2*mu(end);
disp("Cl = " + round(Cl, 3))
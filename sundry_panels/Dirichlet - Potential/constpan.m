clear
close all
clc

%%

alpha = deg2rad(0);
V = 1;

Qinf = [V*cos(alpha); V*sin(alpha)];

N = 300;
coords = genNACA4([2 4 12], 1, N+1);

xe = coords(:,1);
ye = coords(:,2);
xp = xe(2:end) - diff(xe)/2;
yp = ye(2:end) - diff(ye)/2;
thetap = atan2(diff(ye), diff(xe));
np = [-sin(thetap), cos(thetap)];
tp = [cos(thetap), sin(thetap)];
ds = sqrt(diff(xp).^2 + diff(yp).^2);


%%

A = zeros(N+1);
for i = 1:N
    for j = 1:N
        if j ~= i
            A(i, j) = constPhi(xp(i), yp(i), xe(j), ye(j), xe(j+1), ye(j+1), 1);
        else
            A(i, j) = 0.5;
        end
    end
    A(i, N+1) = constPhi(xp(i), yp(i), 1, 0, 1000, 0, 1);
end

A(N+1, :) = [1, zeros(1,N-2), -1, 1];

RHS = -[[xp, yp]*Qinf; 0];

mu = A\RHS;

Qt = diff(mu(1:end-1))./ds;
Cp = 1 - Qt.^2/(Qinf'*Qinf);



figure
hold on
grid on
grid minor
plot(xp(2:end)-diff(xp)/2, Cp)
set(gca, 'ydir', 'reverse')
% xlim([0, 1])
% ylim([-1, 1])

% plot(xp, mu(1:N))

Cl = -2*mu(end);
disp("Cl = " + round(Cl, 3))

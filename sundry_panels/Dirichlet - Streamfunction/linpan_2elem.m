clear
close all
clc

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%%
alpha = deg2rad(0);
Qinf = [cos(alpha); sin(alpha)];

%% element 1
N1 = 101;
wU1 = [0.16 0.16 0.22 0.12];
wL1 = [-0.09 -0.06 -0.12 0.12];
dzTE1 = 0.01;

coords1 = genCST(wU1, wL1, dzTE1, N1)*0.8;
% coords1 = genNACA4([2 2 05], 1, N1);
xp1 = coords1(:,1);
yp1 = coords1(:,2);

thetaTEU1 = -atan(3*sqrt(3)/4*wU1(end)-dzTE1/2);
thetaTEL1 = -atan(3*sqrt(3)/4*wL1(end)+dzTE1/2);
thetaTE1 = 0.5*thetaTEU1 + 0.5*thetaTEL1;

%% element 2
N2 = 101;
wU2 = [0.12 0.12];
wL2 = [-0.12 0.12];
dzTE2 = 0.01;

coords2 = genCST(wU2, wL2, dzTE2, N2)*[cosd(30), -sind(30); sind(30), cosd(30)]*0.2;
% coords2 = genNACA4([2 2 05], 1, N2);
xp2 = coords2(:,1)+0.795;
yp2 = coords2(:,2)-0.025;

thetaTEU2 = -atan(3*sqrt(3)/4*wU2(end)-dzTE2/2);
thetaTEL2 = -atan(3*sqrt(3)/4*wL2(end)+dzTE2/2);
thetaTE2 = 0.5*thetaTEU2 + 0.5*thetaTEL2;

%% combination
xc = [xp1; xp2];
yc = [yp1; yp2];

% left and right panel endpoints
xpl = [xp1(1:N1-1); xp2(1:N2-1)];
xpr = [xp1(2:N1); xp2(2:N2)];
ypl = [yp1(1:N1-1); yp2(1:N2-1)];
ypr = [yp1(2:N1); yp2(2:N2)];

ds = sqrt((xpr-xpl).^2 + (ypr-ypl).^2);
s = [0; cumsum(ds)];

ds12 = [ds(1); 0.5*ds(1:N1-2)+0.5*ds(2:N1-1); ds(N1-1); ds(N1); 0.5*ds(N1:N1+N2-3)+0.5*ds(N1+1:N1+N2-2); ds(N1+N2-2)];

%% AIC
A = zeros(N1+N2+2);

% left and right hats
Psil = zeros(N1+N2, N1+N2-1);
Psir = zeros(N1+N2, N1+N2-1);
% offset for point and panel index mismatch
offsetj = [zeros(1, N1-1), ones(1, N2-1)];
tic
for i = 1:N1+N2
    for j = 1:N1+N2-2
        Psil(i, j+offsetj(j)) = lvs(xc(i), yc(i), xpr(j), ypr(j), xpl(j), ypl(j))*ds(j);
        Psir(i, j+offsetj(j)) = lvs(xc(i), yc(i), xpl(j), ypl(j), xpr(j), ypr(j))*ds(j);
    end
end
toc
% combine left and right influences
A(1:N1+N2, 1:N1+N2-1) = Psil;
A(1:N1+N2, 2:N1+N2) = A(1:N1+N2, 2:N1+N2) + Psir;

% element 1 blunt TE
if dzTE1 ~= 0
    for i = 1:N1
        A(i, 1) = A(i, 1) - 0.5*abs(cos(thetaTE1))*css(xc(i), yc(i), xpl(1), ypl(1), xpr(N1-1), ypr(N1-1))*dzTE1*0.8;
        A(i, 1) = A(i, 1) - 0.5*abs(sin(thetaTE1))*cvs(xc(i), yc(i), xpl(1), ypl(1), xpr(N1-1), ypr(N1-1))*dzTE1*0.8;

        A(i, N1) = A(i, N1) + 0.5*abs(cos(thetaTE1))*css(xc(i), yc(i), xpl(1), ypl(1), xpr(N1-1), ypr(N1-1))*dzTE1*0.8;
        A(i, N1) = A(i, N1) + 0.5*abs(sin(thetaTE1))*cvs(xc(i), yc(i), xpl(1), ypl(1), xpr(N1-1), ypr(N1-1))*dzTE1*0.8;
    end
end
% element 2 blunt TE
if dzTE2 ~= 0
    for i = N1+1:N1+N2
        A(i, N1+1) = A(i, N1+1) - 0.5*abs(cos(thetaTE2))*css(xc(i), yc(i), xpl(N1), ypl(N1), xpr(N1+N2-2), ypr(N1+N2-2))*dzTE2*0.2;
        A(i, N1+1) = A(i, N1+1) - 0.5*abs(sin(thetaTE2))*cvs(xc(i), yc(i), xpl(N1), ypl(N1), xpr(N1+N2-2), ypr(N1+N2-2))*dzTE2*0.2;

        A(i, N1+N2) = A(i, N1+N2) + 0.5*abs(cos(thetaTE2))*css(xc(i), yc(i), xpl(N1), ypl(N1), xpr(N1+N2-2), ypr(N1+N2-2))*dzTE2*0.2;
        A(i, N1+N2) = A(i, N1+N2) + 0.5*abs(sin(thetaTE2))*cvs(xc(i), yc(i), xpl(N1), ypl(N1), xpr(N1+N2-2), ypr(N1+N2-2))*dzTE2*0.2;
    end
end

% element 1 kutta
A(N1+N2+1, :) = [1, zeros(1, N1-2), 1, zeros(1, N2), 0, 0];
% element 2 kutta
A(N1+N2+2, :) = [zeros(1, N1), 1, zeros(1, N2-2), 1, 0, 0];

% element 1 streamfun
A(1:N1, N1+N2+(1:2)) = [-ones(N1, 1), zeros(N1, 1)];
% element 2 streamfun
A(N1+1:N1+N2, N1+N2+(1:2)) = [zeros(N2, 1), -ones(N2, 1)];

PsiA = A(1:N1+N2, 1:N1+N2);

%% RHS
b = [[yc, -xc]*Qinf; 0; 0];

if dzTE1 == 0
    % element 1 extrap
    A(N1, :) = [1, -2, 1, zeros(1, N1-6), -1, 2, -1, zeros(1, N2), 0, 0]; b(N1) = 0;
end
if dzTE2 == 0
    % element 2 extrap
    A(N1+N2, :) = [zeros(1, N1), 1, -2, 1, zeros(1, N2-6), -1, 2, -1, 0, 0]; b(N1+N2) = 0;
end

% figure
% imagesc(A)
% axis equal
% xlim([0, N1+N2+2])
% ylim([0, N1+N2+2])

x = A\b;

g = x(1:N1+N2);
Psi0 = x(N1+N2+1:N1+N2+2);

qt = g;
% qt([1, N]) = -qt([1, N]);
Cp = 1 - qt.^2;
% 
Cl = 2*sum(qt.*[0.5*ds(1); 0.5*ds(1:N1-2)+0.5*ds(2:N1-1); 0.5*ds(N1-1); ...
                0.5*ds(N1); 0.5*ds(N1:N1+N2-3)+0.5*ds(N1+1:N1+N2-2); 0.5*ds(N1+N2-2)]);
disp("Cl = " + Cl)

%% plots

colors = lines(7);

figure
hold on
set(gca, 'ydir', 'reverse')
% plot(xc, Cp, 'linestyle', 'none', 'marker', '.')
plot([xpl(1:N1-1), xpr(1:N1-1)]', [Cp(1:N1-1), Cp(2:N1)]', 'color', colors(1,:))
plot([xpl(N1:N1+N2-2), xpr(N1:N1+N2-2)]', [Cp(N1+1:N1+N2-1), Cp(N1+2:N1+N2)]', 'color', colors(2,:))

% figure('name', 'streamplot', 'units', 'inches', 'papersize', [8 6], 'paperposition', [0 0 8 6])
% hold on
% grid on
% grid minor
% axis equal
% 
% % fimplicit(@(x, y) linPan_PsiPlot_2elem((x-0.5)*cos(-alpha)+y*sin(-alpha)+0.5, (0.5-x)*sin(-alpha)+y*cos(-alpha), ...
% %     xpl, ypl, xpr, ypr, g, Qinf, N1, N2, thetaTE1, dzTE1, thetaTE2, dzTE2)-Psi0(1), [-0.5, 1.5, -1.0, 1.0], 'color', colors(1,:))
% % fimplicit(@(x, y) linPan_PsiPlot_2elem((x-0.5)*cos(-alpha)+y*sin(-alpha)+0.5, (0.5-x)*sin(-alpha)+y*cos(-alpha), ...
% %     xpl, ypl, xpr, ypr, g, Qinf, N1, N2, thetaTE1, dzTE1, thetaTE2, dzTE2)-Psi0(2), [-0.5, 1.5, -1.0, 1.0], 'color', colors(2,:))
% 
% plot(polyshape((xp1-0.5)*cos(alpha)+yp1*sin(alpha)+0.5, (0.5-xp1)*sin(alpha)+yp1*cos(alpha)), ...
%     'facealpha', 1, 'facecolor', 0.5+0.5*colors(1,:))
% plot(polyshape((xp2-0.5)*cos(alpha)+yp2*sin(alpha)+0.5, (0.5-xp2)*sin(alpha)+yp2*cos(alpha)), ...
%     'facealpha', 1, 'facecolor', 0.5+0.5*colors(2,:))
% 
% xlim([-0.5, 1.5])
% xlabel('x/c')
% ylabel('y/c')
% % print('streamplot2.pdf', '-dpdf', '-painters')

%%


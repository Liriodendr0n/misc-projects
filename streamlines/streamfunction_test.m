clear
close all
clc

%% thickness and camber

thickness = 0.20;
camber = 0;
% camber = 3*sqrt(3)/32 * thickness; % zero curvature of lower surface at TE
alpha = deg2rad(10);

ytfun = @(xc) thickness * xc.^0.5 .* (1 - xc) * 3*sqrt(3)/2 + 0.02*xc;
ycfun = @(xc) camber * xc.^1.0 .* (1 - xc) * 4;

npoints = 101;

betas = 2*pi*(0:npoints-1)/(npoints-1) - pi;
xs = (1/2) * (1 - cos(betas));

x = xs;
y = ycfun(xs) + ytfun(xs).*sign(betas)./2;

load foils.mat
x = fliplr(af13.x');
y = fliplr(af13.y');
N = length(x);

coords = [x; y]';
thetas = atan2(diff(y), diff(x));

%% flow solution

[q, g, Vt, Cp, A, r, beta] = panelHS(coords, 1, alpha);

xp = + x.*cos(alpha) + y.*sin(alpha);
yp = - x.*sin(alpha) + y.*cos(alpha);
thetasp = atan2(diff(yp), diff(xp));

%% plot with suction arrows

% figure
% hold on
% axis equal
% plot(xp, yp, 'color', 'black')
% %plot(x(1:end-1) + diff(x)./2, Cp)
% 
% xqp = xp(1:end-1) + diff(xp)./2;
% yqp = yp(1:end-1) + diff(yp)./2;
% 
% kscale = 0.25;
% 
% for i = 1:length(Cp)
%     if Cp(i) <= 0
%         quiver(xqp(i) - kscale*Cp(i).*sin(thetasp(i)).*(1/2 + sign(Cp(i))/2), yqp(i) + kscale*Cp(i).*cos(thetasp(i)).*(1/2 + sign(Cp(i))/2), kscale*Cp(i).*sin(thetasp(i)), -kscale*Cp(i).*cos(thetasp(i)), 0, 'color', [0, 0.4470, 0.7410])
%     else
%         quiver(xqp(i) - kscale*Cp(i).*sin(thetasp(i)).*(1/2 + sign(Cp(i))/2), yqp(i) + kscale*Cp(i).*cos(thetasp(i)).*(1/2 + sign(Cp(i))/2), kscale*Cp(i).*sin(thetasp(i)), -kscale*Cp(i).*cos(thetasp(i)), 0, 'color', [0.8500, 0.3250, 0.0980])
%     end
% end
% 
% figure
% hold on
% plot(x(46:end), -Cp(45:end)+fliplr(Cp(1:44)))
% plot(x(1:end-1), -Cp)
% plot(linspace(0,1,100), 5/pi *sqrt((1-linspace(0,1,100))./linspace(0,1,100)))
% ylim([-1 3])

figure
hold on
axis equal
plot(x(1:end-1)+diff(x)/2, Cp)
set(gca, 'ydir', 'reverse')

%%

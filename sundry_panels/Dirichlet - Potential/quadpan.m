clear
close all
clc

%%

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

xf1 = @(t) interp1(1:N+1, xp, t) + zf(0, -1, t-1).*cos(thetap(1)+pi/2).*ds(1);
yf1 = @(t) interp1(1:N+1, yp, t) + zf(0, -1, t-1).*sin(thetap(1)+pi/2).*ds(1);
  
figure
hold on
axis equal
grid on
grid minor
for i = 1:N
    fplot(@(t) interp1(1:N+1, xp, t) + zf(beta1(i), beta2(i), t-i).*cos(thetap(i)+pi/2).*ds(i), ...
          @(t) interp1(1:N+1, yp, t) + zf(beta1(i), beta2(i), t-i).*sin(thetap(i)+pi/2).*ds(i), [i, i+1], 'linewidth', 1.5)
end
% quiver(xp, yp, cos(np), sin(np), 'color', 'black')
% plot(xp, yp, 'marker', '.', 'color', 'black')


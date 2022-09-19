clear
close all
clc

%% random vectors

u = randn(3,1); v = randn(3,1); w = randn(3,1);
% u = [1; 0; 0]; v = [0; 1; 0]; w = [-1; -1; -0.1];
w = w/norm(w);

%% law of cosines

ca = dot(u,v)/(norm(u)*norm(v));
cb = dot(w,u)/(norm(w)*norm(u));
cc = dot(v,w)/(norm(v)*norm(w));

sa = norm(cross(u,v))/(norm(u)*norm(v));
sb = norm(cross(w,u))/(norm(w)*norm(u));
sc = norm(cross(v,w))/(norm(v)*norm(w));

A = acos((ca-cb*cc)/(sb*sc));
B = acos((cb-ca*cc)/(sa*sc));
C = acos((cc-ca*cb)/(sa*sb));

Oc = A+B+C-pi;

%% simplex formula

Os = 2*atan2(dot(u,cross(v,w)), norm(u)*norm(v)*norm(w) + dot(u,v)*norm(w) + dot(v,w)*norm(u) + dot(w,u)*norm(v));

%% complex formula (works for any n-gon)

% collect normalized vectors
s = [u/norm(u), v/norm(v), w/norm(w)];
% shift
sm1 = circshift(s, 1, 2);
sp1 = circshift(s, -1, 2);

a = dot(sm1,sp1,1);
b = dot(sm1,s,1);
c = dot(s,sp1,1);
d = dot(sm1,cross(s,sp1,1),1);

Oi = 2*pi-sum(atan2(d,b.*c-a));
% correct range
Oi = mod(Oi+2*pi, 4*pi)-2*pi;

% alternate method, only one atan needed, hard to track branches
% Oi2 = -angle(prod(b.*c-a+d*1i));

Oc;
Os
Oi

% p = dot(u,cross(v,w))



% figure
% hold on
% axis equal
% grid on
% grid minor
% view(3)
% % fimplicit3(@(x,y,z) x.^2+y.^2+z.^2-1, [-2, 2], 'FaceAlpha', 0.1, 'edgecolor', 'none')
% quiver3(0,0,0,u(1)/norm(u),u(2)/norm(u),u(3)/norm(u),0)
% quiver3(0,0,0,v(1)/norm(v),v(2)/norm(v),v(3)/norm(v),0)
% quiver3(0,0,0,w(1)/norm(w),w(2)/norm(w),w(3)/norm(w),0)
% plot3([u(1)/norm(u),v(1)/norm(v),w(1)/norm(w),u(1)/norm(u)], [u(2)/norm(u),v(2)/norm(v),w(2)/norm(w),u(2)/norm(u)], [u(3)/norm(u),v(3)/norm(v),w(3)/norm(w),u(3)/norm(u)])
% fimplicit3(@(x,y,z) x.^2+y.^2+z.^2-1, [-2, 2], 'FaceAlpha', 0.1, 'edgecolor', 'none')

bb = 1;

Pz2 = @(aa,u,v) bb * (cos(aa) * u + sin(aa) * v);
aa = fsolve(@(t) Pz2(t, -pi,1)-Pz2(t, pi,1)+2.06, 1.20);

Px = @(u,v) bb*(cos(aa)*sinh(v) .* sin(u) + sin(aa) * cosh(v) .* cos(u));
Py = @(u,v) bb*(-cos(aa)*sinh(v) .* cos(u) + sin(aa) * cosh(v) .* sin(u));
Pz = @(u,v) bb * (cos(aa) * u + sin(aa) * v);

figure
hold on
axis equal
view(3)
fsurf(Px, Py, Pz, [-pi, pi, -1, 1])
% fplot3(@(t) Px(t,1), @(t) Py(t,1), @(t) Pz(t,1), [-pi, pi])

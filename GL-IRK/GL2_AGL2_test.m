clear
close all
clc

%% DISCLAIMER
% looks nice
% doesn't work

% anti-gauss requires higher order information from the true function
% lower order estimates can not produce such information


%% ODE

% y' = lambda*y
lambda = 1;
f = @(t, y) lambda*y;
t0 = 0;
y0 = 1;

%% Method

dt = 0.1;
N = floor(10/dt);

c1 = 1/2;
a1 = 1/2;
b1 = 1;

c2 = [1/2-sqrt(3)/6; 1/2+sqrt(3)/6];
b2 = [1/2, 1/2];
a2 = [1/4, 1/4-sqrt(3)/6; 1/4+sqrt(3)/6, 1/4];

be = [1/2+sqrt(3)/2, 1/2-sqrt(3)/2];
be2 = [1/15, 3/10, 4/15, 3/10, 1/15];

kf = @(k, t, y, a, c) f(t+c*dt, y+a*k*dt) - k;

opts = optimset('Display','off', 'TolFun', 1e-16);

t = [t0, dt*(1:N)];
y = [y0, zeros(1, N)];
% y1 = [y0, zeros(1, N)];
% y2 = [y0, zeros(1, N)];

for i = 1:N
    ks1 = fsolve(@(k) kf(k, t(i), y(i), a1, c1), 1, opts);
    dy1 = dt*b1*ks1;
    
    ks2 = fsolve(@(k) kf(k, t(i), y(i), a2, c2), [1; 1], opts);
    dy2 = dt*b2*ks2;
    
%     y1(i+1) = y1(i) + dy1;
%     y2(i+1) = y2(i) + dy2;
    y(i+1) = y(i) + dy2;
    ke = [y(i); y(i) + a2(1,:)*ks2*dt; y(i) + a1*ks1*dt; y(i) + a2(2,:)*ks2*dt; y(i+1)];
    ye(i+1) = y(i) + dt*be2*ke;
    err(i) = dy2 - dt*be2*ke;
    err2(i) = dy2 - dt*be*ks2;
end

y_exact = y0(1)*exp(lambda*t);

err_global = abs((y_exact(end) - y(end))/y_exact(end))
% err1 = y_exact(end) - y1(end)
% err2 = y_exact(end) - y2(end)

% plot(t, y_exact, t, y)
clear
close all
clc

%% ODE

% y' = lambda*y
lambda = 0.1+1i;
f = @(t, y) lambda*y;
t0 = 0;
y0 = 1;

%% Method

dt = 0.1;
N = floor(10/dt);

c = [0, 2/3];
b = [1/4, 3/4];
a = [1/4, -1/4; 1/4, 5/12];

kf = @(k, t, y) [f(t + c(1)*dt, y + dt*(a(1,1)*k(1) + a(1,2)*k(2))) - k(1); ...
                 f(t + c(2)*dt, y + dt*(a(2,1)*k(1) + a(2,2)*k(2))) - k(2)];

opts = optimset('Display','off', 'TolFun', 1e-16);

t = [t0, dt*(1:N)];
y = [y0, zeros(1, N)];

for i = 1:N
    ks = fsolve(@(k) kf(k, t(i), y(i)), [1; 1], opts);
    y(i+1) = y(i) + dt*b*ks;
end

y_exact = y0(1)*exp(lambda*t);
err = abs(y_exact(end) - y(end));

figure
hold on
plot(t, y_exact)
plot(t, y)


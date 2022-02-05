clear
close all
clc

%% ODE

% y' = lambda*y
lambda = 1i;
f = @(t, y) lambda*y;
t0 = 0;
y0 = 1;
% jacobian
jf = @(t, y) lambda;


%% Method

dt = 0.1;
N = floor(10/dt);

c = [1/2-sqrt(3)/6; 1/2+sqrt(3)/6];
b = [1/2, 1/2];
a = [1/4, 1/4-sqrt(3)/6; 1/4+sqrt(3)/6, 1/4];

be = [1/2+sqrt(3)/2, 1/2-sqrt(3)/2];

kf = @(k, t, y) odeFun(t+c*dt, y+a*k*dt, lambda) - k;

opts = optimset('Display','off', 'TolFun', 1e-16);

t = [t0, dt*(1:N)];
y = [y0, zeros(1, N)];
ye = [y0, zeros(1, N)];

yi = zeros(1, 3*N);
ti = zeros(1, 3*N);

tic
ksold = y0*ones(size(c));
for i = 1:N
    
%     [~, j1] = odeFun(t+c(1)*dt, x+(a(1,1)*k(1)+a12*k2)*dt);
%     [~, j2] = odeFun(x+(a21*k1+a22*k2)*dt);
    
    ks = fsolve(@(k) kf(k, t(i), y(i)), [ksold], opts);
    ksold = ks;
    y(i+1) = y(i) + dt*b*ks;
    
    yi(3*i-2) = y(i);
    yi(3*i-1) = y(i) + a(1,:)*ks*dt;
    yi(3*i-0) = y(i) + a(2,:)*ks*dt;
    
    ti(3*i-2) = t(i);
    ti(3*i-1) = t(i) + c(1)*dt;
    ti(3*i-0) = t(i) + c(2)*dt;
end
toc

tic
opts2 = odeset('abstol', 1e-16);
soln45 = ode45(f, [0 10], y0, opts2);
ti45 = linspace(0, 10, 1000);
yi45 = deval(soln45, ti45);
toc

yiGL = interp1(ti, yi, ti45, 'spline');


y_exact = y0(1)*exp(lambda*ti);

err = abs((y_exact(end) - y(end))/y_exact(end));

figure
hold on
% plot(ti, y_exact)
% plot(t, y-ye)
% plot(t, y)
% plot(ti, yi)
% plot(t, ye)
plot(ti45, yi45)
plot(ti45, yiGL)

function [f, J] = odeFun(t, y, lambda)

% y' = lambda*y

f = lambda*y;

% jacobian
J = lambda;
end


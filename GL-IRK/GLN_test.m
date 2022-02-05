clear
close all
clc

%% derivation

N = 3;

[c, b] = lgwt(N, 0, 1); b = b';

C = zeros(N); V = zeros(N);
for i = 1:N
    for j = 1:N
        V(i, j) = c(i)^(j-1);
        C(i, j) = 1/j * c(i)^j;
    end
end
e = [1./(1:N-1), 0];

a = C/V;
be = e/V;

% main rule order 2N for gauss, 2N-1 for radau, 2N-2 for lobatto
% embedded error estimator order N-1

%% ODE

% y' = lambda*y
lambda = 1i;
f = @(t, y) lambda*y;
t0 = 0;
y0 = 1;
j = @(t, y) lambda+(0*y);

%% Solution

dt = 1;

nt = floor(10/dt);

kf = @(k, t, y) k - f(t+c*dt, y+a*k*dt);

opts = optimset('Display','off', 'TolFun', 1e-16);

t = [t0, dt*(1:nt)];
y = [y0, zeros(1, nt)];

ti = zeros(1, (N+1)*nt);
yi = zeros(1, (N+1)*nt);

ksold = y0*ones(size(c));
tic
for i = 1:nt
%     solve ode system

    %% built in solver (SLOW)
%     ks = fsolve(@(k) kf(k, t(i), y(i)), [ksold], opts);
%     ksold = ks;
    
    %% newton solver (not so slow)
    
    % euler step
    ke = f(t(i), y(i));
    for kj = 1:N
        ks(kj,1) = f(t(i)+c(kj)*dt, y(i) + c(kj)*dt*ke);
    end
    
    Ner = kf(ks, t(i), y(i));
    iters = 0;
%     J = zeros(N);
    while norm(Ner) > 1e-15
        iters = iters+1;
        ji = j(t(i)+c*dt, y(i) + a*ks*dt);
        for Ji = 1:N
            J(Ji,:) = dt*kron(a(Ji,:), ji(Ji));
        end
        J = eye(N) - J;
        k_next = ks - J\Ner;
        ks = k_next;
        Ner = kf(ks, t(i), y(i));
%         iters
    end
    
    
    % take major step
    y(i+1) = y(i) + dt*b*ks;
    
    % step error estimate
    Ser = dt*(b-be)*ks;
    
    % interpolate minor steps
    ti((N+1)*i-N) = t(i);
    yi((N+1)*i-N) = y(i);
    
    ti((N+1)*i-(N-1:-1:0)) = t(i) + c*dt;
    yi((N+1)*i-(N-1:-1:0)) = y(i) + a*ks*dt;
end
ti(end) = t(end);
yi(end) = y(end);
toc

y_exact = y0(1)*exp(lambda*t);
global_err = abs((y_exact(end) - y(end))/y_exact(end))

yifun = @(t) interp1(ti, real(yi), t, 'spline');

figure
hold on
plot(t, real(y))
plot(ti, real(yi))
fplot(yifun, [t(1), t(end)])

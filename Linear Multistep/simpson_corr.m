clear
close all
clc

%%

h = 0.01;
t0 = 0;
t1 = 1;

% lambda = -3.08737802538;
lambda = -1;

y0 = 1;

% y' = lambda*y
f = @(t, y) lambda*y;

%%

t = t0:h:t1;

Y = zeros(size(t));
Y(1) = y0;

% start with ode45
opts = odeset('RelTol', 1e-13);
destart = ode45(f, [0 h], y0, opts);
Y(2) = deval(destart, h);

for n = 2:length(t)-1
    % order 4
    Y(:,n+1) = -4*Y(:,n) + 5*Y(:,n-1) + 2*h*(2*f(t(n), Y(:,n)) + f(t(n-1), Y(:,n-1)));
    Y(:,n+1) = Y(:,n-1) + h/3*(f(t(n+1), Y(:,n+1)) + 4*f(t(n), Y(:,n)) + f(t(n-1), Y(:,n-1)));

%     % order 6 (UNSTABLE)
%     Y(:,n+1) = -18*Y(:,n) + 9*Y(:,n-1) + 10*Y(:,n-2) + 3*h*(3*f(t(n),Y(:,n)) + 6*f(t(n-1),Y(:,n-1)) + f(t(n-2),Y(:,n-2)));
%     Y(:,n+1) = -27/11*Y(:,n) + 27/11*Y(:,n-1) + Y(:,n-2) + 3*h/11*(f(t(n+1),Y(:,n+1)) + 9*f(t(n),Y(:,n)) + 9*f(t(n-1),Y(:,n-1)) + f(t(n-2),Y(:,n-2)));

%     % order 8 (UNSTABLE)
%     Y(:,n+1) = -128/3*Y(:,n) - 36*Y(:,n-1) + 64*Y(:,n-2) + 47/3*Y(:,n-3) + 4*h*(4*f(t(n),Y(:,n)) + 18*f(t(n-1),Y(:,n-1)) + 12*f(t(n-2),Y(:,n-2)) + f(t(n-3),Y(:,n-3)));
%     Y(:,n+1) = -32/5*Y(:,n) + 32/5*Y(:,n-2) + Y(:,n-3) + 6*h/25*(f(t(n+1),Y(:,n+1)) + 16*f(t(n),Y(:,n)) + 36*f(t(n-1),Y(:,n-1)) + 16*f(t(n-2),Y(:,n-2)) + f(t(n-3),Y(:,n-3)));
end

YE = y0*exp(lambda*(t1-t0));
log10(abs((Y(:,end)-YE)/YE))

figure
hold on
fplot(@(t) exp(lambda*t)*y0, [t0, t1])
plot(t, Y, 'marker', '.')


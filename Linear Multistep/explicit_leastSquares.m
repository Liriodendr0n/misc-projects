function [t, y] = explicit_leastSquares(f, tspan, h, N, P, y0)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% start
if N >= 1
    tstart = -(0:N)*h;
%     opts = odeset('reltol', 1e-12);
    start45 = ode45(@(t, y) f(-t, -y), [0 N*h], y0);
    ystart = deval(start45, -tstart);
else
    ystart = y0;
end

%% calculate coefficients b
tic
s = sym(-N:0);

% W = diag(ones(size(s)));
W = diag(1-(s/(N+1)).^2);

J = s'.^(0:P);

a = J'*W*J\J'*W;

for i = 1:N+1
    b(N+2-i) = sum(polyint(flip(a(:,i)')));
end

b = double(b);
toc
%% build t and y histories from start
thist = tspan(1) - h*(N:-1:0);
yhist = ystart;
fhist = zeros(size(ystart));

for i = 1:size(ystart, 2)
    fhist(:,i) = f(thist(i), yhist(:,i));
end
fhist = circshift(fhist, -1, 2);
% fhist = flip(fhist);
% y value at t+0h
yt(:,1) = yhist(:,1);

%% main loop
ts = tspan(1):h:tspan(2);
tic
for ti = 1:length(ts)
    % shift fhist to the right
    fhist = circshift(fhist, 1, 2);
    % set new f(t+0, y(t+0))
    fhist(:,1) = f(ts(ti), yt(:,ti));
    % compute y(t+h) based on history
    ytp1 = yt(:,ti) + h*sum(b.*fhist, 2);
    % update and save y and f and proceed to the next time step
    yt(:,ti+1) = ytp1;
    ft(:,ti+1) = fhist(:,1);
end
toc
%% gather results
t = ts;
y = yt(:,1:end-1);
f = ft;

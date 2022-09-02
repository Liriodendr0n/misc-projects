function [t, y, f] = bpmc2(f, tspan, h, stencil, y0)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% start
if max(stencil) >= 1
    tstart = -(0:max(stencil))*h;

    opts = odeset('reltol', 1e-12);
    start45 = ode45(@(t, y) f(-t, y), [0 (max(stencil)*h)], y0, opts);
    ystart = deval(start45, -tstart);
else
    ystart = y0;
end

%% calculate coefficients b
lagrpol = @(x, x0, j) prod((x - x0(1:end ~= j)')./(x0(j) - x0(1:end ~= j)'), 1);

bi = zeros(1,length(stencil));
gi = zeros(1,length(stencil));
b = zeros(1, max(stencil));

for i = 1:length(stencil)
    if stencil == 0
        bi(i) = 1;
        gi(i) = 1;
        b(stencil(i)+1) = bi(i);
        g(stencil(i)+1) = gi(i);
    else
        bi(i) = integral(@(x) lagrpol(x, stencil, i), -1, 0, 'AbsTol', 1e-15);
        b(stencil(i)+1) = bi(i);
        gi(i) = integral(@(x) lagrpol(x, stencil, i), 0, 1, 'AbsTol', 1e-15);
        g(stencil(i)+1) = gi(i);
    end
end

%% build t and y histories from start
thist = tspan(1) - h*(1:max(stencil)+1);
yhist = ystart;
fhist = zeros(size(ystart));

for i = 1:size(ystart, 2)
    fhist(:,i) = f(thist(i), yhist(:,i));
end
% y value at t+0h
yt(:,1) = yhist(:,1);
%% main loop
ts = tspan(1):h:tspan(2);

for ti = 1:length(ts)
    % shift fhist to the right
    fhist = circshift(fhist, 1, 2);
    % set new f(t+0, y(t+0))
    fhist(:,1) = f(ts(ti), yt(:,ti));
    % compute prediction y'(t+h) based on history
    ytp1p(:,1) = yt(:,ti) + h*sum(b.*fhist, 2);
    % correct prediction twice
    ytp1c = yt(:,ti) + h*sum(g(2:end).*fhist(:,1:end-1), 2) + h*g(1)*f(ts(ti)+h, ytp1p(:,1));
    ytp1 = yt(:,ti) + h*sum(g(2:end).*fhist(:,1:end-1), 2) + h*g(1)*f(ts(ti)+h, ytp1c(:,1));
    % update and save y and f and proceed to the next time step
    yt(:,ti+1) = ytp1;
    ft(:,ti+1) = fhist(:,1);
end
%% gather results
t = ts;
y = yt(:,1:end-1);
f = ft;

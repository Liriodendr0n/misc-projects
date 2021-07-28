clear
close all
clc


c = 1;
N = 101;
coords = genNACA4([2 4 12], c, N);



%% geometry discretization
% surface panel edges
X = coords(:,1);
Y = coords(:,2);

N = length(X)-1;

% control points
Xc = X(1:N) + diff(X)/2;
Yc = Y(1:N) + diff(Y)/2;

% panel lengths
l = sqrt(diff(X).^2 + diff(Y).^2);
sth = diff(Y)./l;
cth = diff(X)./l;

% normal and tangential vectors
n = [-sth; cth];
t = [cth; sth];

% panel global angle
theta = atan2(diff(Y), diff(X));

%% coordinate generation
a = c/2.1;
Xe = Xc - ones(size(Xc))*c/2;
Ye = Yc;

opts = optimset('Diagnostics','off', 'Display','off');
for i = 1:N
    [xsoln] = fsolve(@(X) ellip2cart([Xe(i),Ye(i)],X,a), [1,1], opts);
    mu(i) = xsoln(1);
    nu(i) = xsoln(2);
end

dispN = 10;
muDisp = [0, linspace(0.01,2,dispN-1)];
for i = 1:N
    for j = 1:dispN
        muD(i,j) = mu(i) + muDisp(j)*sign(mu(i));
        nuD(i,j) = nu(i);
    end
end
muD = reshape(muD, 1, N*dispN);
nuD = reshape(nuD, 1, N*dispN);

figure
hold on
scatter(a.*cosh(muD).*cos(nuD), a.*sinh(muD).*sin(nuD))
plot(Xe, Ye)


%% functions
function F = ellip2cart(x,y,a)
F(1) = a*cosh(y(1))*cos(y(2)) - x(1);
F(2) = a*sinh(y(1))*sin(y(2)) - x(2);
end  
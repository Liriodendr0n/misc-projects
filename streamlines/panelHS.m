function [q, g, Vt, Cp, A, r, beta] = panelHS(coords, Vinf, alpha)
%panelHS Hess-Smith sources+vortex panel method
%   Takes airfoil surface coordinates, flight speed, angle of attack
%   Returns source strengths q, vortex strengths g ...
%       surface tangential velocity, surface pressure coefficient

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

% theta(1) = theta(1)-1.5;
% theta(end) = theta(end)-1.5;

%% make r(i,j) and beta(i,j)

r = zeros(N, N+1);
nu = zeros(N, N+1);
beta = zeros(N, N);

for i = 1:N
    for j = 1:N+1
        r(i,j) = sqrt((Xc(i) - X(j))^2 + (Yc(i) - Y(j))^2);
        nu(i,j) = atan2((Yc(i) - Y(j)),(Xc(i) - X(j)));
    end
    beta(i,:) = nu(i,2:end) - nu(i,1:end-1);
    beta(i,i) = pi;
    for j = 1:N
        if beta(i,j) > pi
            beta(i,j) = beta(i,j) - 2*pi;
        end
        if beta(i,j) < -pi
            beta(i,j) = beta(i,j) + 2*pi;
        end
    end
end

%% make A(i,j)
% A(1...N, 1...N)
% tic
A = zeros(N+1, N+1);
%A = gpuArray(A);

for i = 1:N
    for j = 1:N
        A(i,j) = 1/(2*pi) * (sin(theta(i)-theta(j)) * log(r(i,j+1)/r(i,j)) ...
            + cos(theta(i) - theta(j)) * beta(i,j));  
    end
    A(i,N+1) = 1/(2*pi) * ( sum( cos(theta(i)-theta(:)) .* log(r(i,2:end)./r(i,1:end-1))' ...
        - sin(theta(i)-theta(:)).*beta(i,:)'));
end

% A(N+1, 1...N)
for j = 1:N
    A1(j) = sin(theta(1)-theta(j)).*beta(1,j) ...
        - cos(theta(1)-theta(j)).*log(r(1,j+1)./r(1,j));
    AN(j) = sin(theta(N)-theta(j)).*beta(N,j) ...
        - cos(theta(N)-theta(j)).*log(r(N,j+1)./r(N,j));
    A(N+1,j) = 1/(2*pi) * (A1(j) + AN(j));
end

% A(N+1, N+1)
for j = 1:N
    A1Np1(j) = sin(theta(1)-theta(j)).*log(r(1,j+1)./r(1,j)) ...
        + cos(theta(1)-theta(j)).*beta(1,j);
    ANNp1(j) = sin(theta(N)-theta(j)).*log(r(N,j+1)./r(N,j)) ...
        + cos(theta(N)-theta(j)).*beta(N,j);
end
A(N+1,N+1) = 1/(2*pi) * (sum(A1Np1) + sum(ANNp1));
% toc
%% make b(i,j)
for i = 1:N
    b(i,1) = Vinf*sin(theta(i) - alpha);
end
b(N+1,1) = -Vinf*cos(theta(1) - alpha) - Vinf*(cos(theta(N) - alpha));

%% solve Ax = b
% tic
% A = gpuArray(A);
% b = gpuArray(b);

x = A\b;
% x = gather(x);
% toc

%% process solution
% source strengths
q = x(1:N);
% vortex strength
g = x(N+1);

% tangential velocity
for i = 1:N
    Vtf(i) = Vinf*cos(theta(i) - alpha);
    for j = 1:N
        Vtq(i,j) = q(j)/(2*pi) * (sin(theta(i)-theta(j))*beta(i,j) ...
            - cos(theta(i)-theta(j))*log(r(i,j+1)/r(i,j)));
        Vtg(i,j) = g/(2*pi) * (sin(theta(i)-theta(j))*log(r(i,j+1)/r(i,j)) ...
            + cos(theta(i)-theta(j))*beta(i,j));
    end
    Vt(i) = Vtf(i) + sum(Vtq(i,:)) + sum(Vtg(i,:));
end

% pressure coefficient
Cp = 1 - (Vt.^2)./(Vinf^2);

end


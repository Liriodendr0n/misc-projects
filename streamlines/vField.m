function [U, V, beta, r] = vField(coords, plotcoords, q, g, Vinf, alpha)
%VFIELD Summary of this function goes here
%   Detailed explanation goes here

%% geometry discretization
% surface panel edges
X = coords(:,1);
Y = coords(:,2);

Xp = plotcoords(:,1);
Yp = plotcoords(:,2);

N = length(X)-1;


% panel lengths
l = sqrt(diff(X).^2 + diff(Y).^2);
sth = diff(Y)./l;
cth = diff(X)./l;

% normal and tangential vectors
n = [-sth; cth];
t = [cth; sth];

% panel global angle
theta = atan2(gradient(Y,1:length(Y)), gradient(X,1:length(X)));

%% make r(i,j) and beta(i,j)
for i = 1:length(Xp)
    for j = 1:N+1
        r(i,j) = sqrt((Xp(i) - X(j))^2 + (Yp(i) - Y(j))^2);
        nu(i,j) = atan2((Yp(i) - Y(j)),(Xp(i) - X(j))) - 0;
    end
    beta(i,:) = nu(i,2:end) - nu(i,1:end-1);
    for j = 1:N
        if abs(beta(i,j)) > 1.0*pi
            beta(i,j) = beta(i,j) - 2*pi*sign(beta(i,j));
        end
    end
end

%% process solution

% u component
for i = 1:length(Xp)
    Uf(i) = Vinf*cos(0 - alpha);
    for j = 1:N
        Uq(i,j) = q(j)/(2*pi) * (sin(0-theta(j))*beta(i,j) ...
            - cos(0-theta(j))*log(r(i,j+1)/r(i,j)));
        Ug(i,j) = g/(2*pi) * (sin(0-theta(j))*log(r(i,j+1)/r(i,j)) ...
            + cos(0-theta(j))*beta(i,j));
    end
    U(i) = Uf(i) + sum(Uq(i,:)) + sum(Ug(i,:));
end

% v component
for i = 1:length(Yp)
    Vf(i) = Vinf*cos(pi/2 - alpha);
    for j = 1:N
        Vq(i,j) = q(j)/(2*pi) * (sin(pi/2-theta(j))*beta(i,j) ...
            - cos(pi/2-theta(j))*log(r(i,j+1)/r(i,j)));
        Vg(i,j) = g/(2*pi) * (sin(pi/2-theta(j))*log(r(i,j+1)/r(i,j)) ...
            + cos(pi/2-theta(j))*beta(i,j));
    end
    V(i) = Vf(i) + sum(Vq(i,:)) + sum(Vg(i,:));
end

end


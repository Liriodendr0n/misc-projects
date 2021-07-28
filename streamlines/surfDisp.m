function [dispCoords] = surfDisp(coords, thickness)
%surfDisp Displaces surface coordinates outward by a given thickness
%   Detailed explanation goes here

%% geometry discretization
% surface panel edges
X = coords(:,1);
Y = coords(:,2);

N = length(X)-1;

% panel lengths
l = sqrt(diff(X).^2 + diff(Y).^2);
sth = diff(Y)./l;
cth = diff(X)./l;

% normal and tangential vectors
n = [-sth; cth];
t = [cth; sth];

% panel global angle
theta = atan2(diff(Y), diff(X));

%% displace the coordinates

% all points except the last one
Xd(1:N) = X(1:N) + thickness(1:N).*cos(theta(1:N)+pi/2);
Yd(1:N) = Y(1:N) + thickness(1:N).*sin(theta(1:N)+pi/2);

thickness(N+1) = thickness(N) + (thickness(N)-thickness(N-1))*(l(end)/l(end-1));

% final point
Xd(N+1) = X(N+1) + thickness(N+1)*cos(theta(N)+pi/2);
Yd(N+1) = Y(N+1) + thickness(N+1)*sin(theta(N)+pi/2);

dispCoords = [Xd', Yd'];

end


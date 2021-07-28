function [uvw, n] = BiotSavart(xyz, xyzA, xyzB, Gamma, Rc)
%BiotSavart Uses the biot-savart law to compute "induced velocity"
%   A finite-core vortex option exists too

%% distance from point p to vortex filament A --> B

p = xyz;
a = xyzA;
b = xyzB;

n = (b - a)./vecnorm(b - a,2,1);

r = (a - p) - dot(a - p, n, 1).*n;
R = vecnorm(r,2,1);

%% angles from point p to endpoints A and B

pa = a - p;
pb = b - p;

ab = b - a;

alpha = atan2(vecnorm(cross(pa, ab, 1),2,1), dot(pa, ab, 1));
beta = atan2(vecnorm(cross(pb, ab, 1),2,1), dot(pb, ab, 1));

%% induced velocity magnitude and direction

% if R <= 1e-6
%     Vhat = zeros(size(ab));
%     V = 0;
% else
%     Vhat = cross(r, ab)/norm(cross(r, ab));
%     V = -Gamma/(4*pi) * (cos(alpha) - cos(beta)) * R/(R^2 + Rc^2);
% end

Vhat = cross(r, ab, 1)./vecnorm(cross(r, ab, 1), 2, 1);
V = -Gamma/(4*pi) * (cos(alpha) - cos(beta)) .* R./(R.^2 + Rc.^2);
    
uvw = V.*Vhat;

%% error checking
%cross(r, ab)

end


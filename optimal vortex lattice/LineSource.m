function [uvw] = LineSource(r, rA, rB, sigma, Rc)
%LineSource Computes the velocity induced by a line source
%   A finite-core vortex option exists too (WIP)
%%

a = r - rA;
b = r - rB;

s = a - b;
p = a - s.*dot(a, s, 1)./dot(s, s, 1);

uvw = sigma./(4*pi*(dot(p, p, 1) + Rc^2)) .* (p .* (dot(a, s, 1)./(vecnorm(a, 2, 1).*vecnorm(s, 2, 1)) - dot(b, s, 1)./(vecnorm(b, 2, 1).*vecnorm(s, 2, 1))) - ...
      vecnorm(p, 2, 1) .* s ./ vecnorm(s, 2, 1) .* (dot(a, p, 1)./(vecnorm(a, 2, 1).*vecnorm(p, 2, 1)) - dot(b, p, 1)./(vecnorm(b, 2, 1).*vecnorm(p, 2, 1))) ) ;

end
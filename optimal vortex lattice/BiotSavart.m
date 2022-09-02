function [uvw] = BiotSavart(r, rA, rB, Gamma, Rc)
%BiotSavart Uses the biot-savart law to compute "induced velocity"
%   A finite-core vortex option exists too
%%

a = r - rA;
b = r - rB;

uvw = Gamma./(4*pi) * (1./vecnorm(a, 2, 1) + 1./vecnorm(b, 2, 1)) .* cross(a, b, 1)./(vecnorm(a, 2, 1).*vecnorm(b, 2, 1) + dot(a, b, 1) + Rc);
% uvw = Gamma./(4*pi) * (1./vecnorm(a, 2, 1) + 1./vecnorm(b, 2, 1)) .* (a+b)./(vecnorm(a, 2, 1).*vecnorm(b, 2, 1) + dot(a, b, 1) + Rc);



end


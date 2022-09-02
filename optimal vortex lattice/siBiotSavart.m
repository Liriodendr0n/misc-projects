function [uvw] = siBiotSavart(r, rA, t, Gamma, Rc)
%BiotSavart Uses the biot-savart law to compute "induced velocity"
%   A finite-core vortex option exists too
%%

a = r - rA;

uvw = Gamma./(4*pi) * (1./vecnorm(a, 2, 1)) .* cross(a, t, 1)./(vecnorm(a, 2, 1) - dot(a, t, 1) + Rc);


end


function [uvw] = iBiotSavart(r, rA, t, Gamma, Rc)
%BiotSavart Uses the biot-savart law to compute "induced velocity"
%   A finite-core vortex option exists too
%%

a = r - rA;

uvw = Gamma./(2*pi) * cross(a, t, 1)./(dot(a, a, 1) - dot(a, t, 1).^2 + Rc);


end


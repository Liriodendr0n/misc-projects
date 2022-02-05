function [V] = linesource(z, z0, z1, a0, a1)
%cubicsource Computes the induced velocity around a cubic panel with cubic strength
%   Detailed explanation goes here

%% panel polynomial coeffs

c0 = z0;
c1 = z1 - z0;

th = angle(z1-z0);

%% induced velocity

C = (a1*(z-c0) + a0*c1) * log((c1+c0-z)/(c0 - z));

V = exp(2i*th)*conj(-(C + a1*c1)/(2*pi*abs(c1)));

end




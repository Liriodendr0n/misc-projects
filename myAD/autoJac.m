function [J] = autoJac(f, x)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

xx = DS;
Nvar = length(x);
for i = 1:length(x)
    xx(i) = DSvar(Nvar, 1, i, x(i));
end

J = DSjac(f(xx));

end


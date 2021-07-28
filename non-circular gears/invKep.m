function [E] = invKep(M, e)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

options = optimoptions('fsolve', 'display', 'off');

fun = @(E) E - e.*sin(E) - M;

E = fsolve(fun, M,options);

end


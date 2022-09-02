function [dy] = lagrPolyDer(j, X, x)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

k = length(X)-1;

sumterm = zeros(k+1, length(x));

for m = 0:k
    if m ~= j
        sumterm(m+1,:) = lagrPoly(j, X, x)./(x - X(m+1));
    else
        sumterm(m+1,:) = zeros(size(x));
    end
end

dy = sum(sumterm, 1);

end
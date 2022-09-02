function [y] = lagrPoly(j, X, x)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

k = length(X)-1;

prodterm = ones(k+1, length(x));

for m = 0:k
    if m ~= j
        prodterm(m+1,:) = (x - X(m+1))./(X(j+1) - X(m+1));
    else
        prodterm(m+1,:) = ones(size(x));
    end
end

y = prod(prodterm, 1);

% y = prod((x - X(1:end ~= j+1)')./(X(j+1) - X(1:end ~= j+1)'), 1);

end
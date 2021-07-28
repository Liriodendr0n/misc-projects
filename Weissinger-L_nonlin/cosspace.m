function [Xk] = cosspace(x1, x2, n)
%chebspace spaces n points akin to chebyshev nodes on (x1, x2)
%   Detailed explanation goes here

k = 1:n;

xk = cos((2*k - 2)/(2*n-2) * pi);

Xk = x2/2*(1 - xk) + x1/2*(1 + xk);


end


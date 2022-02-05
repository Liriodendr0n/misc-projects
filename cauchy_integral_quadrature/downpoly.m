function [p] = downpoly(N, x)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


nCk = @(n, k) factorial(n)./(factorial(k).*factorial(n-k));

m = sym(0:N);

pN = nCk(2*floor(m/2), floor(m/2))./(4.^floor(m/2) .* (-1).^(m+1));

p = dot(pN, x.^(N:-1:0));

end


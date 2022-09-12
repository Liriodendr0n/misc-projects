function [A, dA] = interpMatrix3(xv, xq, alpha, beta)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

N = length(xv);
M = length(xq);
Q = repmat(xq(:), [1,length(xv)])';
R = repmat(xv(:), [1,length(xv)])';

W = (1-Q).^alpha .* (1+Q).^beta;
dW = beta*(1-Q).^alpha.*(1+Q).^(beta-1) - alpha*(1-Q).^(alpha-1).*(1+Q).^beta;

R(logical([triu(ones(N,N-1)), zeros(N,1)])) = R(logical(triu(ones(N), 1)));
R = R(:,1:end-1);
% rows of R are roots of lagrange polynomials

% weights to satisfy kronecker delta property
w = prod(1./(xv(:)-R), 2);

% 3d indexing :)
xq3 = reshape(xq, [1, 1, M]);
l = reshape(prod((xq3-R), 2), [N, M]);
dl = reshape(sum(1./(xq3-R), 2), [N, M]);

L = (l.*w)';
dL = (l.*dl.*w)';

A = L.*W';
dA = L.*dW' + dL.*W';

end


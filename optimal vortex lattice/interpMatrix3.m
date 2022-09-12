function [A, dA] = interpMatrix3(xv, xq, alpha, beta)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

xq(xq>1) = 1;
xv(xv>1) = 1;
xq(xq<-1) = -1;
xv(xv<-1) = -1;

N = length(xv);
M = length(xq);
Q = repmat(xq(:), [1,length(xv)])';
R = repmat(xv(:), [1,length(xv)])';

% jacobi weight function
W = (1-Q).^alpha .* (1+Q).^beta;

if alpha == 0 && beta ~= 0
    dW = beta.*(1+Q).^(beta-1);
elseif alpha ~= 0 && beta == 0
    dW = -alpha.*(1-Q).^(alpha-1);
elseif alpha == 0 && beta == 0
    dW = zeros(N,M);
else
    dW = beta*(1-Q).^alpha.*(1+Q).^(beta-1) - alpha*(1-Q).^(alpha-1).*(1+Q).^beta;
end

% eliminating diagonal from matrix of roots
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


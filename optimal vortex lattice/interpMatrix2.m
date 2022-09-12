function [A, dA] = interpMatrix2(xv, xq, alpha, beta)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

N = length(xv);
Q = repmat(xq(:), [1,length(xv)])';
R = repmat(xv(:), [1,length(xv)])';

W = (1-Q).^alpha .* (1+Q).^beta;
dW = beta*(1-Q).^alpha.*(1+Q).^(beta-1) - alpha*(1-Q).^(alpha-1).*(1+Q).^beta;

R(logical([triu(ones(N,N-1)), zeros(N,1)])) = R(logical(triu(ones(N), 1)));
R = R(:,1:end-1);
% rows of R are roots of lagrange polynomials

% weights to satisfy kronecker delta property
w = prod(1./(xv(:)-R), 2);

% interpolating polynomial coefficients
P = cellfun(@poly, num2cell(R, 2), 'UniformOutput', false);
P = cellfun(@times, P, num2cell(w), 'UniformOutput', false);
dP = cellfun(@polyder, P, 'UniformOutput', false);

% interpolation matrices
A = cell2mat(cellfun(@polyval, P, num2cell(Q, 2), 'UniformOutput', false));
dA = cell2mat(cellfun(@polyval, dP, num2cell(Q, 2), 'UniformOutput', false));

% weighted interpolation matrices
dA = (A.*dW + dA.*W)';
A = (A.*W)';

end


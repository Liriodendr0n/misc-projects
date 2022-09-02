function [RHS] = vorlatRHS(alpha, n)
%weisslRHS Computes the right hand side vector RHS
%   Includes nonlinear section data

vinf = [cos(alpha), 0, sin(alpha)]';

N = size(n, 2);

% loop over control points
% for j = 1:N
%     % flow penetration is allowed at the control points
%     % this accounts for the lift curve not being 2*pi*alpha
%     RHS(j,1) = Vtransp(j) - dot(vinf, n(:,j));
% end

RHS1 = ones(N,1);
J = find(RHS1);

RHS = dot(repmat(vinf, 1, N), n(:,J), 1)';

end


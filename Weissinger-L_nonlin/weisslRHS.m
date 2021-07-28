function [RHS] = weisslRHS(alpha, n, Vtransp)
%weisslRHS Computes the right hand side vector RHS
%   Includes nonlinear section data

vinf = [cos(alpha), 0, sin(alpha)]';

N = length(Vtransp);

% loop over control points
% for j = 1:N
%     % flow penetration is allowed at the control points
%     % this accounts for the lift curve not being 2*pi*alpha
%     RHS(j,1) = Vtransp(j) - dot(vinf, n(:,j));
% end

RHS1 = ones(N,1);
J = find(RHS1);

RHS = Vtransp(J) - dot(repmat(vinf, 1, N), n(:,J), 1)';

end


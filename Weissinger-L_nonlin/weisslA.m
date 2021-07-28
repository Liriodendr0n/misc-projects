function [A, a1, a2, a3] = weisslA(xyzC, xyzV, n)
%weisslA Computes the influence matrix A
%   Planar wings only, but planform can be whatever

N = size(xyzC, 2);

% loop over all control points
% for i = 1:N
%     
%     % loop over all vortices
%     for j = 1:N 
%         % build horseshoe vortices
%         dTV = [1e9; 0; 0];
%         
%         % left trailing vortex (infty to left endpoint)
%         a1 = BiotSavart(xyzC(:,i), xyzV(:,j)+dTV, xyzV(:,j), 1, 0);
%         % bound vortex (left to right)
%         a2 = BiotSavart(xyzC(:,i), xyzV(:,j), xyzV(:,j+1), 1, 0);
%         % right trailing vortex (right endpoint to infty)
%         a3 = BiotSavart(xyzC(:,i), xyzV(:,j+1), xyzV(:,j+1)+dTV, 1, 0);
%         
%         % sum the influences of each line vortex to form the horseshoe
%         A(i, j) = dot((a1+a2+a3), n(:,i));
%         
%     end
% end

dTV = [1e9; 0; 0];

[I, J] = find(ones(N));

a1 = BiotSavart(xyzC(:,I), xyzV(:,J)+dTV, xyzV(:,J), 1, 0);
a2 = BiotSavart(xyzC(:,I), xyzV(:,J), xyzV(:,J+1), 1, 0);
a3 = BiotSavart(xyzC(:,I), xyzV(:,J+1), xyzV(:,J+1)+dTV, 1, 0);
A1 = dot((a1+a2+a3), n(:,I), 1);

A = reshape(A1, N, N);

end


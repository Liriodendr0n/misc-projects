function [Atfz] = vorlatAtrefftz(xyzD, xyzV, n)
%weisslA Computes the influence matrix A
%   Planar wings only, but planform can be whatever

N = size(xyzD, 2);

% for i = 1:N
%     % loop over all horseshoe vortices
%     for j = 1:N 
%         % build horseshoe vortices
%         dTV = [1e6; 0; 0];
%         
%         % left trailing vortex (infty to left endpoint)
%         a1 = BiotSavart(xyzD(:,i), xyzV(:,j)+dTV, xyzV(:,j), 1, 0);
%         % bound vortex (left to right)
%         a2 = BiotSavart(xyzD(:,i), xyzV(:,j), xyzV(:,j+1), 1, 1e-6);
%         % right trailing vortex (right endpoint to infty)
%         a3 = BiotSavart(xyzD(:,i), xyzV(:,j+1), xyzV(:,j+1)+dTV, 1, 0);
%         
%         % sum the influences of each line vortex to form the horseshoe
%         AindV(i, j) = dot((a1+a2+a3), n(:,i));
%         
%     end
% end

dTV = [1e9; 0; 0];

[I, J] = find(ones(N));

a1 = BiotSavart(xyzD(:,I), xyzV(:,J)+dTV, xyzV(:,J)-dTV, 1, 0);
a2 = BiotSavart(xyzD(:,I), xyzV(:,J+1)-dTV, xyzV(:,J+1)+dTV, 1, 0);

AindV1 = dot((a1+a2), n(:,I), 1);

Atfz = reshape(AindV1, N, N);

end


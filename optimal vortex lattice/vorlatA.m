function [A] = vorlatA(xyzC, xyzVL, xyzVR, n)
%weisslA Computes the influence matrix A
%   Planar wings only, but planform can be whatever

N = size(xyzC, 2);

tTV = repmat([1; 0; 0], 1, N);

% collect index pairs for vectorized evaluation
[I, J] = find(ones(N));

a1 = siBiotSavart(xyzC(:,I), xyzVL(:,J), tTV(:,J), 1, 1e-14);
a2 = BiotSavart(xyzC(:,I), xyzVL(:,J), xyzVR(:,J), 1, 1e-14);
a3 = siBiotSavart(xyzC(:,I), xyzVR(:,J), tTV(:,J), -1, 1e-14);

A1 = dot((a1+a2+a3), n(:,I), 1);

A = reshape(A1, N, N);

end


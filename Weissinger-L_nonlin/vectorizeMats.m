clear
close all
clc

%%

alpha = deg2rad(5);
% ||v_infty|| is assumed to be 1

xyz = [0; 0; 0];
N = 1000;

AR = 3;
S = 1;
lambda = 0.2;
Lambda = deg2rad(33);
Gamma = deg2rad(0);


b = sqrt(AR*S);
c = S/b;
%% build geometry
[xyzV, xyzC, xyzD, cCs, cVs, n] = wingGeom(xyz, N, AR, S, lambda, Lambda, Gamma);

tic
dTV = [1e9; 0; 0];

[I, J] = find(ones(N));

[a11, n1] = BiotSavart(xyzC(:,I), xyzV(:,J)+dTV, xyzV(:,J), 1, 0);
a21 = BiotSavart(xyzC(:,I), xyzV(:,J), xyzV(:,J+1), 1, 0);
a31 = BiotSavart(xyzC(:,I), xyzV(:,J+1), xyzV(:,J+1)+dTV, 1, 0);
A11 = dot((a11+a21+a31), n(:,I), 1);

A1 = reshape(A11, N, N);
toc

tic
[A, a1, a2, a3] = weisslA(xyzC, xyzV, n);
toc

RHS = Vtransp(J) - dot(vinf, n(:,J), 1);


function [v, c] = chord_LEKutta(N, Gamma)
%spanTip Summary of this function goes here
%   Detailed explanation goes here

%%

% polynomial: a = 0, b = 0
% chordwise: a = 0.5, b = -0.5
% spanwise: a = 0.5, b = 0.5

a = 1;
b = -0.5;

%% Recurrence Coefficients

ab = a + b;
ii = (2:N-1)';
abi = 2*ii + ab;
aa = [(b - a)/(2 + ab)
      (b^2 - a^2)./((abi - 2).*abi)
      (b^2 - a^2)./((2*N - 2+ab).*(2*N+ab))];
bb = [2*sqrt( (1 + a)*(1 + b)/(ab + 3))/(ab + 2) ; 
      2*sqrt(ii.*(ii + a).*(ii + b).*(ii + ab)./(abi.^2 - 1))./abi];

% weight function integral
IW = 2^(a+b+1)*((gamma(a+1)*gamma(b+1))/gamma(a+b+2));


% Gamma = 0;

bb(N-1) = sqrt(1+Gamma)*bb(N-1);


%% Jacobi Matrix

J = diag(bb, 1) + diag(aa) + diag(bb, -1);

% special case
if N == 1
    J = aa(1);
end

%% Radau modifications

Rp = -1;

bbN = 2*sqrt(N.*(N + a).*(N + b).*(N + ab)./((2*N + ab).^2 - 1))./(2*N + ab);
JR = J-(Rp*eye(N));
eN = [zeros(N-1, 1); 1];

delta = JR\([bb; bbN].^2.*eN);
dN = delta(N);

aaNp1 = Rp + dN;

J = [J, bbN*eN; bbN*eN', aaNp1];


%% Golub-Welsch

[V, D] = eig(J);

v = diag(D);
w = IW*V(1,:).^2;

%% Cauchy Nodes

W = @(t) (1-t).^a .* (1+t).^b;

% precomputed principal value of weight function
C = @(t) pi;

Q = @(t) w*(1./(t-v));

f = @(t) C(t) - Q(t);

c = zeros(N, 1);

% figure
% hold on
% fplot(f, [-1, 1])
% ylim([-1, 1])
% xlim([-2, 2])

vint = [v; 1.01];

for ci = 1:N+1
    c(ci, 1) = fzero(f, [vint(ci)+1e-6, vint(ci+1)-1e-6]);
end

% figure
% hold on
% fplot(W, [-1, 1])
% fplot(f, [-1, 1])
% ylim([-2, 2])

%% Interval Change

v = flip((1+v')/2);
c = flip((1+c')/2);

end


clear
close all
clc


N = 2;


%%

% polynomial: a = 0, b = 0
% chordwise: a = 0.5, b = -0.5
% spanwise: a = 0.5, b = 0.5

a = 0;
b = 0;

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

% Gamma = 1+1/(N+2);
Gamma = 0
bb(N-1) = sqrt(1+Gamma)*bb(N-1);
%% Jacobi Matrix

J = diag(bb, 1) + diag(aa) + diag(bb, -1);

% special case
if N == 1
    J = aa(1);
end

%% Golub-Welsch

[V, D] = eig(J);

v = diag(D);
w = IW*V(1,:).^2;

%% Cauchy

W = @(t) (1-t).^a .* (1+t).^b;

% precomputed principal value of weight function
C = @(t) log((1+t)./(1-t));

Q = @(t) w*(1./(t-v));

f = @(t) C(t) - Q(t);

c = zeros(N, 1);

vint = [-1; v; 1];

for ci = 1:N+1
    c(ci, 1) = fzero(f, [vint(ci)+1e-12, vint(ci+1)-1e-12]);
end

figure
hold on
fplot(f, [-1, 1])
fplot(W, [-1, 1])
fplot(@(t) 0.*t, [-1, 1])
ylim([-2, 2])

c
clear
% close all
clc

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

syms x w theta

%% calculate coefficients b

N = 16;
P = 8;

tic
s = sym(-N:0);

% W = diag(ones(size(s)));
W = diag(1-(s/(N+1)).^2);

J = s'.^(0:P);

a = J'*W*J\J'*W;

for i = 1:N+1
    b(i) = sum(polyint(flip(a(:,i)')));
end

toc

alpha = [zeros(1, N), -1, 1];

beta = [b 0];
% beta = [flip([1.99846268500208 -1.28551740667502 -0.217557223320753 0.797815926658251 -0.293203981664557]), 0]
k = length(beta)-1;

h(w) = sum(alpha.*w.^(0:k))./sum(beta.*w.^(0:k));

s(theta) = exp(1i*theta);

f = matlabFunction(h(s));

xt = @(t) real(f(t));
yt = @(t) imag(f(t));

% figure
hold on
grid on
grid minor
axis equal
color = lines(7);
plot(xt(linspace(-pi, pi, 3000)), yt(linspace(-pi, pi, 3000)))

% pretty(flip(beta(1:end-1)));

%%


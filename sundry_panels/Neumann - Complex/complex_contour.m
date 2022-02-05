clear
close all
clc

%%

syms z0 z1 th0 th1 t

syms c0 c1 c2 c3

c = [c0 c1 c2 c3];

%% cubic curve

gamma(t) = c0 + c1*t + c2*t^2 + c3*t^3;

eqn = [gamma(0) == z0; gamma(1) == z1; ...
       subs(diff(gamma, t), t, 0) == abs(z1-z0)*exp(1i*th0); subs(diff(gamma, t), t, 1) == abs(z1-z0)*exp(1i*th1)];

soln = solve(eqn, [c]);

c0 = soln.c0;
c1 = soln.c1;
c2 = soln.c2;
c3 = soln.c3;

% gamma(t) = c0 + c1*t + c2*t^2 + c3*t^3;

%% integrand

syms z zeta f0 f1

% green's function
G(z, zeta) = 1/sym(2*pi) * 1/(z-zeta);


% value
syms a0 a1 a2 a3
f(t) = a0 + a1*t;

dg(t) = diff(gamma, t);

ifun = simplify(G(z, gamma)*abs(dg(t))*f(t))




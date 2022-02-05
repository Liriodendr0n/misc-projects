clear
close all
clc

%%

N = 1000;

coords = genNACA4([2 4 12], 1, N);
[q, g, Vt, Cp, A, b, r, beta] = panelHS(coords, 1, 0);

tic
A\b;
toc


testfun = @(x) A*x - b;

tic
[Xsoln] = fsolve(testfun, 0.0*[q; g]);
toc

mean(abs(Xsoln - [q; g]))
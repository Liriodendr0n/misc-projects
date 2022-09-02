clear
close all
clc

%%

syms h z

alpha = [1, 0, -1];
beta = [1/3, 4/3, 1/3];

rho = alpha*z.^(0:2).';
sigma = beta*z.^(0:2).';
clear
close all
clc

%%

%% system
% x = DSvar(3, 2, 1, 4);
% y = DSvar(3, 2, 2, 5);
% z = DSvar(3, 2, 3, 14);
% 
% u = [x, y; 1, z];
% 
% F = sqrt(x+y)./sqrt(z-y);
% tic
% % F = exp(x.*sin(y))./(z.*tan(x)-cos(y));
% 
% Vf = DSval(F);
% Gf = DSgrad(F);
% Hf = DShess(F);
% toc

%% inverse

x = DSvar(1, 2, 1, 1);

A = [randn, x; randn, randn];

I = inv(A);

Afd = A.f;
Afd_dx = Afd + [0, 1e-9; 0, 0];

(inv(Afd_dx) - inv(Afd))./1e-9;
-inv(Afd)*[0, 1; 0, 0]*inv(Afd);
I.df.f

2*inv(Afd)*[0, 1; 0, 0]*inv(Afd)*[0, 1; 0, 0]*inv(Afd);
I.df.df;

X = ainit(1, 2);
Y = inv([Afd(1,1), X; Afd(2,1), Afd(2,2)]);
dY = adiff(Y, 1);
ddY = adiff(Y, 2);

[Y(1,1).c(2), Y(1,2).c(2); Y(2,1).c(2), Y(2,2).c(2)]
[Y(1,1).c(3), Y(1,2).c(3); Y(2,1).c(3), Y(2,2).c(3)];

%% matmul

% A = [1, 2; 3, 4];
% B = [4, 5; 6, 7];
% 
% x = DSvar(1, 1, 1, 0);
% 
% AD = [A(1,1), A(1,2)+x; A(2,1), A(2,2)];
% BD = [B(1,1), B(1,2)+x; B(2,1), B(2,2)];
% 
% CD = AD*BD;
% 
% dx = 1e-9;
% dX = [0, dx; 0, 0];
% 
% ((A+[0, dx; 0, 0])*(B+[0, dx; 0, 0]) - A*B)/dx
% CD.df


%% eigen

% A = randn(2);
% x = DSvar(1, 2, 1, A(1,2));
% 
% AD = [A(1,1), x; A(2,1), A(2,2)];
% 
% dx = 1e-9;
% Adx = A + [0 dx; 0 0];
% 
% [u, l, v] = eig(A);
% [udx, ldx, vdx] = eig(Adx);
% 
% (ldx-l)/dx;
% 
% lD = eig(AD);
% lD.df
% 
% trace((u(:,1)*v(:,1)'/(v(:,1)'*u(:,1)))*[0 1; 0 0])
% trace((u(:,2)*v(:,2)'/(v(:,2)'*u(:,2)))*[0 1; 0 0])


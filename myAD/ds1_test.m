clear
close all
clc

%%

%% system

u = DSvar(2, 3, 1, 0);
v = DSvar(2, 3, 2, 0);

x = DSvar(3, 2, 1, 4);
y = DSvar(3, 2, 2, 5);
z = DSvar(3, 2, 3, 14);

% x = DSvar(3, 2, 1, 4+u);
% y = DSvar(3, 2, 2, 5+v);
% z = DSvar(3, 2, 3, 14);

% u = [x, y; 1, z];
tic
F = sqrt(x+y)./sqrt(z-y);
toc
% tic
% F = exp(x.*sin(y))./(z.*tan(x)-cos(y));


Vf = DSval(F);
Gf = DSgrad(F);
Hf = DShess(F);

Hf0 = [-0.00308641975308642 0 -0.00308641975308642;0 0.0123456790123457 -0.0123456790123457;-0.00308641975308642 -0.0123456790123457 0.00925925925925926];
Hf1 = [-0.00308641971879287 6.85871056130272e-11 -0.00308641978737997;6.85871056130272e-11 0.0123456794238683 -0.0123456793552812;-0.00308641978737997 -0.0123456793552812 0.00925925956790124];
dHf = (Hf1-Hf0)/1e-7;

partial(Hf, [1 0])
% partial(Hf, 0)
% Hf.f
% % toc

%% nested

% x = DSvar(2, 1, 1, 3);
% y = DSvar(2, 1, 2, 4);

% c = DSvar(1, 1, 1, 3);
% x = DSvar(2, 1, 1, c);
% y = DSvar(2, 1, 2, c+1);
% 
% z = sin(x).*sqrt(y);
% 
% partial(z, [0, 0])
% partial(z, [1, 0])
% partial(z, [0, 1])

%% inverse

% x = DSvar(1, 2, 1, 1);
% A = DSconst(1, 2, randn(10));
% 
% A(4,7) = x;
% 
% I = inv(A);
% 
% partial(I, 1);

%% matmul

% x0 = 1;
% 
% A = [1, 2; 3, 4];
% B = [4, 5; 6, 7];
% 
% x = DSvar(1, 1, 1, 0);
% 
% AD = [A(1,1)+x, A(1,2)+x; A(2,1), A(2,2)];
% BD = [B(1,1), B(1,2)+x; B(2,1), B(2,2)];
% 
% CD = AD*BD;
% 
% dx = 1e-9;
% dX = [0, dx; 0, 0];
% 
% ((A+[dx, dx; 0, 0])*(B+[0, dx; 0, 0]) - A*B)/dx
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


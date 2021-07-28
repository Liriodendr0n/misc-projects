clear
close all
clc

%%

alpha = rad2deg(1/(2*pi));

c = 2;

N = 1000;
k = 0:N;

%% pistolesi

X = -cos(pi*(0:N)/N);
xVp = X(1:end-1) + 1/4*(X(2:end) - X(1:end-1));
xCp = X(1:end-1) + 3/4*(X(2:end) - X(1:end-1));

% build A
for i = 1:N
    for j = 1:N
        Ap(i, j) = 1/(2*pi*(xCp(i) - xVp(j)));
    end
end

% build b
for j = 1:N
    bp(j, 1) = sind(alpha);
end

% solve
Gp = Ap\bp;

errP = abs((Gp./(xCp' - xVp'))./(2/pi * sqrt((1-xVp')./(1+xVp'))) - 1);

Lp = Gp./(xCp' - xVp');
Mp = sum(Gp .*(xVp' + 0.5));

%% Stark?

xVs = -cos((2*(1:N)-1)/(2*N+1) *pi);
xCs = -cos((2*(1:N)-0)/(2*N+1) *pi);

Ws = pi/(2*N+1) * sin((2*(1:N)-1)/(2*N+1) *pi);

% build A
for i = 1:N
    for j = 1:N
        As(i, j) = 1/(2*pi*(xCs(i) - xVs(j)));
    end
end

% build b
for j = 1:N
    bs(j, 1) = sind(alpha);
end

% solve
Gs = As\bs;

errS = abs((Gs./(xCs' - xVs'))./(2/pi * sqrt((1-xVs')./(1+xVs'))) - 1);

Ls = Gs./(xCs' - xVs');
Ms = sum(Gs .*(xVs' + 0.5));

%% plot


sum(Ws'.*Ls./(xCs'-xVs'))


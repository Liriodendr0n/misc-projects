clear
close all
clc

%% Saddle-Free Newton Optimization

a = 1;
b = 100;
f = @(x, y) (a-x).^2 + b.*(y-x.^2).^2

x0 = -3;
y0 = -2;

X = DSvar(2, 2, 1, x0);
Y = DSvar(2, 2, 2, y0);

XY = [X; Y];

xx = x0;
yy = y0;
figure
hold on
axis equal
fcontour(@(x, y) log10(f(x, y)), [-4.5, 4.5], 'linecolor', 'black', 'meshdensity', 300)

tic
for i = 1:3
    
    F = f(X, Y);
    G = DSgrad(F);
    H = DShess(F);
    
    [V, D] = eig(H);
    Hplus = V*abs(D)*inv(V);
    
    XY = XY - H\G;
    X = XY(1);
    Y = XY(2);
    xx = [xx, DSval(X)];
    yy = [yy, DSval(Y)];
    plot(xx, yy, 'marker', '.')
    pause(0.1)
    DSval(F)
end
toc


[DSval(X), DSval(Y)];
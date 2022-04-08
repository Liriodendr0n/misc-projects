clear
close all
clc

load foils.mat

splitind = ceil(size(af17,1)/2);
xs = linspace(0, 1, 100);

yL = interp1(af17.x(splitind:end), af17.y(splitind:end), xs);
yU = interp1(af17.x(1:splitind), af17.y(1:splitind), xs);

yC = (yU+yL)/2;

cMaxOld = max(abs(yC))*sign(max(yC));

cMaxNew = 0.01;

yLnew = 2*yC*(cMaxNew./cMaxOld) - yU;

figure
hold on
axis equal
plot(xs, yL, xs, yU, xs, yLnew)

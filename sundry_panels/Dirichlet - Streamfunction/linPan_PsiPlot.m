function [Psi] = linPan_PsiPlot(x, y, xp, yp, g, Qinf, thetaTE, dzTE)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

N = length(xp);
M = length(x);

ds = sqrt(diff(xp).^2 + diff(yp).^2);

x = x';
y = y';


% xy = ones(M, N);
% xy(:,[1,N]) = zeros(M, 2);
% 
% [I1, J1] = find(xy);
% [I2, J2] = find(1-xy);

A = zeros(M, N);
for i = 1:M
    for j = 2:N-1
        Psi1 = lvs(x(i), y(i), xp(j-1), yp(j-1), xp(j), yp(j))*ds(j-1);
        Psi2 = lvs(x(i), y(i), xp(j+1), yp(j+1), xp(j), yp(j))*ds(j);
        A(i, j) = Psi1 + Psi2;
    end
end

for i = 1:M
    A(i, 1) = lvs(x(i), y(i), xp(2), yp(2), xp(1), yp(1))*ds(1);
    A(i, N) = lvs(x(i), y(i), xp(N-1), yp(N-1), xp(N), yp(N))*ds(N-1);
end

% trailing edge gap
if dzTE ~= 0
    for i = 1:M
        A(i, 1) = A(i, 1) - 0.5*abs(cos(thetaTE))*css(x(i), y(i), xp(1), yp(1), xp(N), yp(N))*dzTE;
        A(i, 1) = A(i, 1) - 0.5*abs(sin(thetaTE))*cvs(x(i), y(i), xp(1), yp(1), xp(N), yp(N))*dzTE;

        A(i, N) = A(i, N) + 0.5*abs(cos(thetaTE))*css(x(i), y(i), xp(1), yp(1), xp(N), yp(N))*dzTE;
        A(i, N) = A(i, N) + 0.5*abs(sin(thetaTE))*cvs(x(i), y(i), xp(1), yp(1), xp(N), yp(N))*dzTE;
    end
end

Psi = sum(A'.*g, 1) + sum([-y, x]'.*Qinf, 1);

end


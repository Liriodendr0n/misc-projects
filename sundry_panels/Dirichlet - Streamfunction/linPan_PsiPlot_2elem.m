function [Psi] = linPan_PsiPlot_2elem(x, y, xpl, ypl, xpr, ypr, g, Qinf, N1, N2, thetaTE1, dzTE1, thetaTE2, dzTE2)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

M = length(x);

A = zeros(M, N1+N2);

ds = sqrt((xpr-xpl).^2 + (ypr-ypl).^2);

x = x';
y = y';

% left and right hats
Psil = zeros(M, N1+N2-1);
Psir = zeros(M, N1+N2-1);
% offset for point and panel index mismatch
offsetj = [zeros(1, N1-1), ones(1, N2-1)];
for i = 1:M
    for j = 1:N1+N2-2
        Psil(i, j+offsetj(j)) = lvs(x(i), y(i), xpr(j), ypr(j), xpl(j), ypl(j))*ds(j);
        Psir(i, j+offsetj(j)) = lvs(x(i), y(i), xpl(j), ypl(j), xpr(j), ypr(j))*ds(j);
    end
end

% combine left and right influences
A(1:M, 1:N1+N2-1) = Psil;
A(1:M, 2:N1+N2) = A(1:M, 2:N1+N2) + Psir;

% element 1 blunt TE
if dzTE1 ~= 0
    for i = 1:M
        A(i, 1) = A(i, 1) - 0.5*abs(cos(thetaTE1))*css(x(i), y(i), xpl(1), ypl(1), xpr(N1-1), ypr(N1-1))*dzTE1;
        A(i, 1) = A(i, 1) - 0.5*abs(sin(thetaTE1))*cvs(x(i), y(i), xpl(1), ypl(1), xpr(N1-1), ypr(N1-1))*dzTE1;

        A(i, N1) = A(i, N1) + 0.5*abs(cos(thetaTE1))*css(x(i), y(i), xpl(1), ypl(1), xpr(N1-1), ypr(N1-1))*dzTE1;
        A(i, N1) = A(i, N1) + 0.5*abs(sin(thetaTE1))*cvs(x(i), y(i), xpl(1), ypl(1), xpr(N1-1), ypr(N1-1))*dzTE1;
    end
end
% element 2 blunt TE
if dzTE2 ~= 0
    for i = 1:M
        A(i, N1+1) = A(i, N1+1) - 0.5*abs(cos(thetaTE2))*css(x(i), y(i), xpl(N1), ypl(N1), xpr(N1+N2-2), ypr(N1+N2-2))*dzTE2;
        A(i, N1+1) = A(i, N1+1) - 0.5*abs(sin(thetaTE2))*cvs(x(i), y(i), xpl(N1), ypl(N1), xpr(N1+N2-2), ypr(N1+N2-2))*dzTE2;

        A(i, N1+N2) = A(i, N1+N2) + 0.5*abs(cos(thetaTE2))*css(x(i), y(i), xpl(N1), ypl(N1), xpr(N1+N2-2), ypr(N1+N2-2))*dzTE2;
        A(i, N1+N2) = A(i, N1+N2) + 0.5*abs(sin(thetaTE2))*cvs(x(i), y(i), xpl(N1), ypl(N1), xpr(N1+N2-2), ypr(N1+N2-2))*dzTE2;
    end
end

Psi = sum(A'.*g, 1) + sum([-y, x]'.*Qinf, 1);

end


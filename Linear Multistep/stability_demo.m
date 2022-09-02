clear
close all
clc

syms x1 x2 x3 x4 real

x = [x1; x2; x3; x4];
f(x) = [x3.*x2.^2 - 1./x3.^2; -2./x3 .*x1.*x2; x1; x2];

A(x) = jacobian(f, x);

Anum = matlabFunction(A);
fnum = matlabFunction(f);
eigsnum = matlabFunction(eig(A));

opts = odeset('reltol', 1e-6);
soln45 = ode45(@(t, x) fnum(x(1), x(2), x(3), x(4)), [0 5*pi], [0, 1.2, 1, 0], opts);

t = soln45.x;
x1 = soln45.y(1,:);
x2 = soln45.y(2,:);
x3 = soln45.y(3,:);
x4 = soln45.y(4,:);

eigvals = zeros(4, length(t));
for ti = 1:length(t)
    Ai = Anum(x1(ti), x2(ti), x3(ti), x4(ti));
%     eigvals(:,ti) = roots(charpoly(Ai));
    eigvals(:,ti) = eigsnum(x1(ti), x2(ti), x3(ti), x4(ti));
end

figure
hold on
axis equal
plot(real(eigvals)', imag(eigvals)')
plot(x3.*cos(x4), x3.*sin(x4), 'marker', '.')

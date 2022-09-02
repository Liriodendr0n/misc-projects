clear
close all
clc

syms x1 x2 x3 x4 real

% % r theta orbit IC = (rdot thetadot r theta)
% f = @(t, y) [y(3).*y(2).^2 - 1./y(3).^2; -2./y(3) .*y(1).*y(2); y(1); y(2)];
% y0 = [0; 1.1; 1; 0];

x = [x1; x2; x3; x4];
f = [x(3).*x(2).^2 - 1./x(3).^2; -2./x(3) .*x(1).*x(2); x(1); x(2)];

Apol(x) = jacobian(f, x);



% % x y orbit IC = (vx vy rx ry)
% f = @(t, y) [-y(3)./((y(3).^2 + y(4).^2).^(3/2)); -y(4)./((y(3).^2 + y(4).^2).^(3/2)); y(1); y(2)];
% y0 = [0; 1.4; 1; 0];

x = [x1; x2; x3; x4];
f = [-x(3)./((x(3).^2 + x(4).^2).^(3/2)); -x(4)./((x(3).^2 + x(4).^2).^(3/2)); x(1); x(2)];

Acart(x) = jacobian(f, x);

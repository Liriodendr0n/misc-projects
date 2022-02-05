clear
close all
clc

%%

syms x y t real

psi1(x, y) = -1/sym(2*pi) * atan(x/y);
psi2(x, y) = -1/sym(2*pi) * atan2(x,y);
psi3(x, y) = -1/sym(2*pi) * atan2(y,-x);

theta1 = atan(x./y);
theta2 = atan((x-1)./y);
r1 = x.^2 + y.^2;
r2 = (x-1).^2 + y.^2;

Psi(x, y) = sign(y)*int(psi1(x-t, y), t, 0, 1);
Psi2(x, y) = -((y.*log(r1./r2))/2 + theta2.*(x - 1) - x.*theta1)/(2*pi) + (0.5*(0.5+0.5*sign(y)) - 0.5);

Psi1(x, y) = ((x-1)*atan2(1-x,y)+(y*log((x-1)^2+y^2))/2) - ((x-0)*atan2(0-x,y)+(y*log((x-0)^2+y^2))/2);

% fcontour(Psi(x,y))
fcontour(Psi2)
% fcontour(psi2)
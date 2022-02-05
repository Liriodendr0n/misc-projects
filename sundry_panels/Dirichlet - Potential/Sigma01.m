function [Phi] = Sigma01(x0, y0, x1, y1, x2, y2, sigma0, sigma1)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


panel_ang = -atan2(y2-y1, x2-x1);
panel_len = sqrt((x2-x1).^2 + (y2-y1).^2);

x = cos(panel_ang).*x0 - sin(panel_ang).*y0 - x1.*cos(panel_ang)+y1.*sin(panel_ang);
y = sin(panel_ang).*x0 + cos(panel_ang).*y0 - x1.*sin(panel_ang)-y1.*cos(panel_ang);


x1 = 0;
x2 = panel_len;
theta1 = atan((x-x1)./y);
theta2 = atan((x-x2)./y);
r1 = (x-x1).^2 + y.^2;
r2 = (x-x2).^2 + y.^2;

Phi = sigma0/(4*pi) * ((x-x1)*log(r1) - (x-x2)*log(r2) + 2*y*(theta1 - theta2)) + sigma1/(8*pi) * ((x^2-x1^2-y^2)*log(r1) - (x^2-x2^2-y^2)*log(r2) + 4*x*y*(theta1-theta2) - 2*x*(x2-x1));

end


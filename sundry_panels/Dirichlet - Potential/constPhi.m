function [Phi] = constPhi(x0, y0, x1, y1, x2, y2, mu0)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% translate and rotate

panel_ang = -atan2(y2-y1, x2-x1);
panel_len = sqrt((x2-x1).^2 + (y2-y1).^2);

x = (cos(panel_ang).*x0 - sin(panel_ang).*y0 - x1.*cos(panel_ang)+y1.*sin(panel_ang))/panel_len;
y = (sin(panel_ang).*x0 + cos(panel_ang).*y0 - x1.*sin(panel_ang)-y1.*cos(panel_ang))/panel_len;


x1 = 0;
x2 = 1;
theta1 = atan((x-x1)./y);
theta2 = atan((x-x2)./y);
r1 = (x-x1).^2 + y.^2;
r2 = (x-x2).^2 + y.^2;

%% calculate
Phi = - mu0/(2*pi)*(theta1 - theta2); %- mu1/(4*pi)*(2*x.*(theta1 - theta2) + y.*log(r2./r1)) + mu2/(2*pi) * ((x.^2-y.^2).*(theta2-theta1) - x.*y.*log(r2./r1) + y.*(x1-x2));

end


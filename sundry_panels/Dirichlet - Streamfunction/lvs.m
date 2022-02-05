function [Psi] = lvs(x0, y0, x1, y1, x2, y2)
%lvs linear vortex streamfunction
%   Detailed explanation goes here

%% translate and rotate

panel_ang = -atan2(y2-y1, x2-x1);
panel_len = sqrt((x2-x1).^2 + (y2-y1).^2);

x = (cos(panel_ang).*x0 - sin(panel_ang).*y0 - x1.*cos(panel_ang)+y1.*sin(panel_ang))./panel_len;
y = (sin(panel_ang).*x0 + cos(panel_ang).*y0 - x1.*sin(panel_ang)-y1.*cos(panel_ang))./panel_len;

%% bits

theta1 = atan(x./y);
theta2 = atan((x-1)./y);
r1 = x.^2 + y.^2;
r2 = (x-1).^2 + y.^2;
h = x.^2 - y.^2;

%% calculate

if x == 0 && y == 0
    Psi = 0;
elseif x == 1 && y == 0
    Psi = 1/(4*pi);
else
    Psi = -1/(4*pi) * (h/2.*log(r1) - (h-1)/2.*log(r2) + 2*x.*y.*(theta1-theta2) - x);
end

end


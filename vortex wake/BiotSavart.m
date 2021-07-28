function [u, v] = BiotSavart(x, xp, y, yp, Gp, r0)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


r = sqrt((x-xp).^2 + (y-yp).^2);

if r < r0
    r = r0.^2 ./r;
end

if r ~= 0
    theta = atan2((y-yp), (x-xp));

    V = Gp./(2*pi*r);

    u = V.*cos(theta + pi/2);
    v = V.*sin(theta + pi/2);
else
    u = 0;
    v = 0;
end

end


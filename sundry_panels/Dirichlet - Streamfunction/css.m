function [Psi] = css(x0, y0, x1, y1, x2, y2)
%lvs constant source streamfunction
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

%% calculate

if y == 0
    Psi = 0;
else
    
    % jumps along y = 0 (in plane of panel)
%     Psi = -((y.*log(r1./r2))/2 + theta2.*(x - 1) - x.*theta1)/(2*pi);

    % patched jump?
    Psi = -((y.*log(r1./r2))/2 + theta2.*(x - 1) - x.*theta1)/(2*pi) *sign(y);

    % jumps along x = 0,1 for y < 0 (below plane of panel)
%     Psi = 1/(2*pi) * (((x-1).*atan2(1-x,y)+(y.*log(r2)/2)) - ((x-0)*atan2(0-x,y)+(y*log(r1)/2)));

    % point source
%     Psi = -1/(2*pi*panel_len) * atan2(y, -x);
end


end


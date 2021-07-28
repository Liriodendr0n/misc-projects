clear
clc
close all

I = [1 2 1];


tspan = [0 30];

w0 = [1 0.1 0.2];

opts = odeset('RelTol', 1e-6);
soln = ode45(@(t, w) Euler(t, w, I), tspan, w0, opts);

figure
hold on
plot(soln.x, soln.y(1,:))
plot(soln.x, soln.y(2,:))
plot(soln.x, soln.y(3,:))
%plot(soln.x, sqrt(soln.y(1,:).^2 + soln.y(2,:).^2 + soln.y(3,:).^2))

figure
hold on
view(3)
axis equal
%quiver3(zeros(size(soln.x)), zeros(size(soln.x)), zeros(size(soln.x)), ...
%    soln.y(1,:), soln.y(2,:), soln.y(3,:))
plot3(soln.y(1,:), soln.y(2,:), soln.y(3,:))



function [dwdt] = Euler(t, w, I)

dwdt(1,:) = -(I(3) - I(2))/I(1) .* w(2)*w(3);
dwdt(2,:) = -(I(1) - I(3))/I(2) .* w(3)*w(1);
dwdt(3,:) = -(I(2) - I(1))/I(3) .* w(1)*w(2);

end


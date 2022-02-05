function [coords] = genCST(wL, wU, dz, N)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

C = @(xc) 3*sqrt(3)/2 * (xc).^0.5 .* (1-xc).^1;
B = @(xc, n, k) factorial(n)./(factorial(k).*factorial(n-k)) .* (xc).^(k) .* (1-xc).^(n-k);

beta = 2*pi*(0:N-1)/(N-1) - pi;
x = (1/2) * (1 - cos(beta'));

if mod(N, 2) == 0
    xU = x(1:N/2);
    xL = x(N/2+1:N);
else
    xU = x(1:(N+1)/2-1);
    xL = x((N+1)/2:N);
end

yU = 0.5*sum(B(xU, length(wU)-1, (0:length(wU)-1)).*wU, 2).*C(xU) - 0.5*dz*xU;
yL = 0.5*sum(B(xL, length(wL)-1, (0:length(wL)-1)).*wL, 2).*C(xL) + 0.5*dz*xL;

y = [yU; yL];

coords = [x, y];

end


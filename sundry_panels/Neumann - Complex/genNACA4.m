function [xout, yout] = genNACA4(code, npoints)
%UNTITLED2 Summary of this function goes here

%% geometric parameters
chord = 1;
c = chord;
m = code(1)/100;     % maximum camber [%]
p = code(2)/10;    % position of maximum camber [tenths]
t = code(3)/100;     % maximum thickness [%]

%% cosine space points
beta = 2*pi*(0:npoints-1)/(npoints-1) - pi;
x = (c/2) * (1 - cos(beta));

%% thickness function
ytFun = @(xc) 5*t*c*(0.2969*(xc).^(1/2) - 0.1260*(xc).^1 - 0.3516*(xc).^2 + 0.2843*(xc).^3 - 0.1036*(xc).^4);

%% camber function
ycFun = @(xc) m/(p^2) .* (2*p*(xc) - xc.^2) .* heaviside(p - xc) ...
    + m/(1-p)^2 .* ((1-2*p) + 2*p*xc - xc.^2) .* heaviside(xc - p);

%% camberline slope function
thetaFun = @(xc) atan2((2*m)/(p^2) * (p-xc) .* heaviside(p - xc) ...
    + (2*m)/((1-p)^2) * (p-xc) .* heaviside(xc - p), 1);

yc = ycFun(x/c);
yt = ytFun(x/c);
theta = thetaFun(x/c);

%% create surface
X = x - yt.*sin(theta).*sign(sin(beta));

Y = yc - yt.*cos(theta).*sign(sin(beta));

xout = X';
yout = Y';

coords = [X', Y'];

end


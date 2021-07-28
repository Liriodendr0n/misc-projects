function [xyzV, xyzC, xyzD, cCs, cVs, ns] = wingGeomKuch(xyz, N, AR, S, lambda, Lambda, Gamma)
%wingGeom Defines the geometry of a general trapezoidal wing
%   xyz is the position of the wing apex (center tip of the LE)
%   N the number of horseshoe vortices, control, downwash points
%   AR = aspect ratio, S = wing area, lambda = tip taper ratio
%   Lambda = quarter-chord sweep angle, Gamma = dihedral angle

%	vortex points are located at c/4
%   control points are located at 3c/4
%   downwash points are located at c/2
%   a vector of chords is provided too

%% planform definition
% wingspan
b = sqrt(AR*S);

% root and tip chords (taper)
cr = 2*b/(AR*(lambda + 1));
ct = lambda*cr;
cfun = @(y) cr - 2*(cr-ct)/b *abs(y);

% root and tip x displacements (sweep)
dxr = 0;
dxt = b/2 * tan(-Lambda);

% root and tip z displacements (dihedral)
dzr = 0;
dzt = b/2 * tan(Gamma);

%% kuchemann locus

LambdaK = Lambda/((1 + ((2*pi*cos(Lambda))/(pi*AR))^2)^(1/4));

K = (1 + ((2*pi*cos(Lambda))/(pi*AR))^2)^(pi/(4*(pi + 2*abs(LambdaK))));

lambfun = @(y) (sqrt(1 + (2*pi * tan(LambdaK+1e-4)/(LambdaK+1e-4) * abs(y)./cfun(y)).^2) - ...
    2*pi * tan(LambdaK+1e-4)/(LambdaK+1e-4) * abs(y)./cfun(y)) - ...
    (sqrt(1 + (2*pi * tan(LambdaK+1e-4)/(LambdaK+1e-4) * (b/2 - abs(y))./cfun(y)).^2) - ...
    2*pi * tan(LambdaK+1e-4)/(LambdaK+1e-4) * (b/2 - abs(y))./cfun(y));

f = @(y) cr/4 + abs(y).*tan(Lambda) - cfun(y)/4 .*(1 - (1/K).*(1 + 2.*lambfun(y).*(LambdaK/pi)));

%% geometry definition

etas = cosspace(-1, 1, N+1);

% vortices
dyVs = etas * b/2;
cVs = cfun(dyVs);
dxVs = f(dyVs);
%dxVs = cr/4 + abs(dyVs)*tan(Lambda);
dzVs = tan(Gamma) * abs(dyVs);

% control points
dyCs = (dyVs(1:end-1) + dyVs(2:end))/2;
cCs = cfun(dyCs');
dxCs = f(dyCs) + cfun(dyCs)/2;
%dxCs = cr/4 + abs(dyCs)*tan(Lambda) + cfun(dyCs)/2;
dzCs = tan(Gamma) * abs(dyCs);

% downwash points
dyDs = dyCs;
cDs = cCs;
%dxDs = dxCs - cCs/2;
%dxDs = f(dyCs) + cfun(dyCs)/4;
dxDs = 1000*ones(size(dyDs));
dzDs = tan(Gamma) * abs(dyDs);

% normal vectors
ns = [0*ones(1,N); -sin(Gamma)*sign(dyCs); cos(Gamma)*ones(1,N)];

%% final geometry

xyzV = xyz + [dxVs; dyVs; dzVs];
xyzC = xyz + [dxCs; dyCs; dzCs];
xyzD = xyz + [dxDs; dyDs; dzDs];


end


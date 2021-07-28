function [xyzV, xyzC, xyzD, cCs, cVs, ns] = wingGeom(xyz, N, AR, S, lambda, Lambda, Gamma)
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

% root and tip x displacements (sweep)
dxr = 0;
dxt = b/2 * tan(-Lambda);

% root and tip z displacements (dihedral)
dzr = 0;
dzt = b/2 * tan(Gamma);

%% geometry definition

etas = cosspace(-1, 1, N+1);

% vortices
dyVs = etas * b/2;
cVs = cr - 2*(cr-ct)/b * abs(dyVs);
dxVs = dxr + 2*(dxr-dxt)/b * abs(dyVs) + cr/4;
dzVs = tan(Gamma) * abs(dyVs);

% control points
dyCs = (dyVs(1:end-1) + dyVs(2:end))/2;
cCs = cr - 2*(cr-ct)/b * abs(dyCs');
dxCs = (dxVs(1:end-1) + dxVs(2:end))/2 + cCs'/2;
dzCs = tan(Gamma) * abs(dyCs);

% downwash points (trefftz plane for now)
dyDs = dyCs;
cDs = cr - 2*(cr-ct)/b * abs(dyDs');
%dxDs = dxCs - cDs'/4;
dxDs = 1e4*ones(size(dyDs));
dzDs = tan(Gamma) * abs(dyDs);

% normal vectors
ns = [0*ones(1,N); -sin(Gamma)*sign(dyCs); cos(Gamma)*ones(1,N)];

%% final geometry

xyzV = xyz + [dxVs; dyVs; dzVs];
xyzC = xyz + [dxCs; dyCs; dzCs];
xyzD = xyz + [dxDs; dyDs; dzDs];


end


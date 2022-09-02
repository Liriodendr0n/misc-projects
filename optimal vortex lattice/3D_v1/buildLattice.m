function [xyzL, xyzR, xyzC, xyzTL, xyzTR, ns] = buildLattice(cornerPts, xiV, xiC, etaV, etaC)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

N1 = @(xi, eta) (1-xi).*(1-eta);
N2 = @(xi, eta) (1-xi).*(1+eta);
N3 = @(xi, eta) (1+xi).*(1+eta);
N4 = @(xi, eta) (1+xi).*(1-eta);

NxiC = length(xiC);
NxiV = length(xiV);
NetaC = length(etaC);
NetaV = length(etaV);

xiT = xiV;
etaT = etaV;

% control points
[xiC, etaC] = meshgrid(xiC, etaC);
xiC = reshape(xiC, [1, NxiC*NetaC]);
etaC = reshape(etaC, [1, NxiC*NetaC]);
xyzC = [cornerPts(1,:)*[0.25*N1(xiC, etaC) ; 0.25*N2(xiC, etaC) ; 0.25*N3(xiC, etaC) ; 0.25*N4(xiC, etaC)]; ...
        cornerPts(2,:)*[0.25*N1(xiC, etaC) ; 0.25*N2(xiC, etaC) ; 0.25*N3(xiC, etaC) ; 0.25*N4(xiC, etaC)]; ...
        cornerPts(3,:)*[0.25*N1(xiC, etaC) ; 0.25*N2(xiC, etaC) ; 0.25*N3(xiC, etaC) ; 0.25*N4(xiC, etaC)]];

% bound vortex endpoints
[xiV, etaV] = meshgrid(xiV, etaV);
xiV = reshape(xiV, [1, NxiV*NetaV]);
etaV = reshape(etaV, [1, NxiV*NetaV]);
xyzV = [cornerPts(1,:)*[0.25*N1(xiV, etaV) ; 0.25*N2(xiV, etaV) ; 0.25*N3(xiV, etaV) ; 0.25*N4(xiV, etaV)]; ...
        cornerPts(2,:)*[0.25*N1(xiV, etaV) ; 0.25*N2(xiV, etaV) ; 0.25*N3(xiV, etaV) ; 0.25*N4(xiV, etaV)]; ...
        cornerPts(3,:)*[0.25*N1(xiV, etaV) ; 0.25*N2(xiV, etaV) ; 0.25*N3(xiV, etaV) ; 0.25*N4(xiV, etaV)]];
% split left and right
xyzL = xyzV(:,setdiff(1:NxiV*NetaV, 1:NetaV:NxiV*NetaV));
xyzR = xyzV(:,setdiff(1:NxiV*NetaV, NetaV:NetaV:NxiV*NetaV));

% trailing vortices
[xiT, etaT] = meshgrid(ones(size(xiT)), etaT);
xiT = reshape(xiT, [1, NxiV*NetaV]);
etaT = reshape(etaT, [1, NxiV*NetaV]);
xyzT = [cornerPts(1,:)*[0.25*N1(xiT, etaT) ; 0.25*N2(xiT, etaT) ; 0.25*N3(xiT, etaT) ; 0.25*N4(xiT, etaT)]; ...
         cornerPts(2,:)*[0.25*N1(xiT, etaT) ; 0.25*N2(xiT, etaT) ; 0.25*N3(xiT, etaT) ; 0.25*N4(xiT, etaT)]; ...
         cornerPts(3,:)*[0.25*N1(xiT, etaT) ; 0.25*N2(xiT, etaT) ; 0.25*N3(xiT, etaT) ; 0.25*N4(xiT, etaT)]];
% split left and right
xyzTL = xyzT(:,setdiff(1:NxiV*NetaV, 1:NetaV:NxiV*NetaV));
xyzTR = xyzT(:,setdiff(1:NxiV*NetaV, NetaV:NetaV:NxiV*NetaV));

% placeholder normal vectors (zhat only)
ns = [zeros(1,NxiV*(NetaV-1)); zeros(1,NxiV*(NetaV-1)); ones(1,NxiV*(NetaV-1))];

end


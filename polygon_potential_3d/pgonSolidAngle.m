function [Omega] = pgonSolidAngle(C, p)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% collect normalized vectors

s = (C-p)./vecnorm(C-p, 2, 1);
% shift
sm1 = circshift(s, -1, 2);
sp1 = circshift(s, 1, 2);

%% 10.48550/arXiv.1205.1396

% a = dot(sm1,sp1,1);
% b = dot(sm1,s,1);
% c = dot(s,sp1,1);
% d = dot(sm1,cross(s,sp1,1),1);
% 
% % slow but guaranteed to be accurate
% omicron = 2*pi-sum(atan2(d,b.*c-a));
% 
% % fast but not guaranteed, branch issues with complex argument
% % omicron = 2*pi-angle(prod(b.*c-a+d*1i));

%% tangent vectors

% taum = sm1 - s;
% taup = sp1 - s;
% omicron = 2*pi - sum(atan2(vecnorm(cross(taum, taup, 1),2,1),dot(taum, taup, 1)));

%% sum of triangles

s0 = [1; 0; 0].*ones(size(s));
omicron = sum(2*atan2(dot(s,cross(sp1,s0,1)), 1 + dot(s,sp1,1) + dot(sp1,s0,1) + dot(s0,s,1)));

% correct range
Omega = mod(omicron+2*pi, 4*pi)-2*pi;
% Omega = omicron;

end


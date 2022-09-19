function [Omega] = curveSolidAngle(C, p)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% collect normalized vectors

s = (C-p)./vecnorm(C-p, 2, 1);
% shift
sm1 = circshift(s, -1, 2);
sp1 = circshift(s, 1, 2);

a = dot(sm1,sp1,1);
b = dot(sm1,s,1);
c = dot(s,sp1,1);
d = dot(sm1,cross(s,sp1,1),1);

omciron = 2*pi-sum(atan2(d,b.*c-a));
% correct range
Omega = mod(omciron+2*pi, 4*pi)-2*pi;

end


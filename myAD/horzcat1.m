function [out] = horzcat1(a, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    out = a;
if nargin == 1
    return
elseif nargin > 2
    out = horzcat1(out, horzcat1(varargin{1}, varargin{2:end}));
else
    out = horzcat(out, varargin{:});
end
    

end


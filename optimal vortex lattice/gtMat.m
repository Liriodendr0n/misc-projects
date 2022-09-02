function [A] = gtMat(xv, t)
%gtMat does stuff
%   

A = zeros(length(xq), length(xv));


for i = 1:length(t)
    J(:,i) = jacobiP(i-1, 0.5, -0.5, 2*xiV-1);
end



for i = 1:length(xq)
    for j = 1:length(xv)
        B(i,j) = lagrPoly(j-1, xv, xq(i));     
    end
end




end
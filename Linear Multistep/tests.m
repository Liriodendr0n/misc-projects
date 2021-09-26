clear
clc
close all

x0 = [0 2 5 7];

lagrpol = @(x, x0, j) prod((x - x0(1:end ~= j)')./(x0(j) - x0(1:end ~= j)'), 1);

bi = zeros(1,length(x0));
b = zeros(1, max(x0));
for i = 1:length(x0)
    bi(i) = integral(@(x) lagrpol(x, x0, i), -1, 0);
    b(x0(i)+1) = bi(i);
end

b
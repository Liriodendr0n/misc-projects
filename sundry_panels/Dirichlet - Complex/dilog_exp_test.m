clear
close all
clc

%%

f = @(z) polylog(2, 1 - exp(z));
g = @(z) floor(imag(z/(2*pi))-1/2)+1;



Li2 = @(z) polylog(2, z) + (2i*pi).*log(z) .* g(z);


figure
hold on
view(2)
fsurf(@(x, y) real(Li2(1 - exp(x+1i*y))), 'edgecolor', 'none')
% fsurf(@(x, y) real(exp(x+1i*y)), 'edgecolor', 'none')
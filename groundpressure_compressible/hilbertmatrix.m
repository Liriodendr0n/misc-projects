clear
close all
clc

Ns = 1:30;

for i = 1:length(Ns)
    
    N = Ns(i);
    A = hilb(N);
    check = A\A;

    err(i) = abs(det(check)-1);
    
end


plot(Ns, err)
set(gca, 'Yscale', 'log')
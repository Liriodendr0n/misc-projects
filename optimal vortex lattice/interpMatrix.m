function [A, dA, IA] = interpMatrix(xv, xq, alpha, beta)
%interpMatrix Creates a "Lagrange Transition Matrix"
%   It interpolates points xq from the lagrange basis on xv

A = zeros(length(xq), length(xv));
dA = zeros(length(xq), length(xv));
IA = zeros(length(xq), length(xv));

for i = 1:length(xq)
    for j = 1:length(xv)
        A(i,j) = lagrPoly(j-1, xv, xq(i))*(1-xq(i))^alpha*(1+xq(i))^beta;
        dA(i,j) = lagrPolyDer(j-1, xv, xq(i))*(1-xq(i))^alpha*(1+xq(i))^beta + ...
                  lagrPoly(j-1, xv, xq(i))*( beta*(1-xq(i))^alpha*(1+xq(i))^(beta-1) - alpha*(1-xq(i))^(alpha-1)*(1+xq(i))^beta);
%         IA(i,j) = integral(@(t) lagrPoly(j-1, xv, t).*(1-t).^alpha*(1+t).^beta, -0.999, xq(i), 'arrayValued', true, 'abstol', 1e-12, 'reltol', 1e-12);
    end
end

end
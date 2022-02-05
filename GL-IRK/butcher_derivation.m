clear
close all
clc

%%

% nodes
c = [-1/sym(3)*sqrt(5+2*sqrt(sym(10)/7)), -1/sym(3)*sqrt(5-2*sqrt(sym(10)/7)), 0, 1/sym(3)*sqrt(5-2*sqrt(sym(10)/7)), 1/sym(3)*sqrt(5+2*sqrt(sym(10)/7))];
% c = (1+[-sym(sqrt(3))/3, sym(sqrt(3))/3])/2;

% weights
b = [(322-13*sym(sqrt(70)))/900, (322+13*sym(sqrt(70)))/900, 128/225, (322+13*sym(sqrt(70)))/900, (322-13*sym(sqrt(70)))/900];


for i = 1:5
    for j = 1:5
        V(i, j) = c(i)^(j-1);
        C(i, j) = 1/j * c(i)^j;
        N(i, j) = 1/i;
    end
end
D = diag(b);
e = [sym(1), sym(1)./(2:5-1), 0];

%% Gaussian
A = C*inv(V);
be = simplify(e/V);


%% Radau IA
% A = inv(D)*inv(V')*(N-C)'*D;
% be = simplify(e/V);


%% Radau IIA
% A = C*inv(V);
% be = simplify(e/V);


%% Lobatto IIIA
% A = C*inv(V);
% be = simplify(e/V);


%% Lobatto IIIB
% A = inv(D)*inv(V')*(N-C)'*D;
% be = simplify(e/V);


%% Lobatto IIIC
% for i = 1:3
%     for j = 1:2
%         C2(i, j) = (c(i)^j)/j;
%     end
%     C2(i, 1) = C2(i, 1) - b(1);
% end
% for i = 1:2
%     for j = 1:2
%         V2(i, j) = c(i+1)^(j-1);
%     end
% end
% A(:, 1) = b(1)*ones(3,1);
% A(:,2:3) = C2*inv(V2);
% be = simplify(e/V);


%%


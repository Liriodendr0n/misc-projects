clear
close all
clc

%%

A = rand(3, 2);
B = rand(2, 5);

A*B;

% a = reshape(A, [numel(A),1]);
% b = reshape(B.', [numel(B),1]);


C = zeros(size(A,1), size(B,2));
tic
for i = 1:size(A, 1)
    for j = 1:size(B, 2)
        C(i,j) = 0;
        for k = 1:size(A, 2)
            C(i,j) = C(i,j) + A(i,k).*B(k,j);
        end
    end
end
toc
C;
%{
    ForwardSub(L, b) solves L * x = b, 
    where L is a lower triangular matrix. (Forward Substitution)
%}

function [x] = ForwardSub(L, b)

% Find size of L
row = length(L);

% Initialize solution
x = zeros(row, 1);

x(1) = b(1) / L(1, 1);

% Forward substitution
for i = 2 : row
    x(i) = b(i);
    for j = 1 : i - 1
        x(i) = x(i) - L(i, j) * x(j);
    end
    x(i) = x(i) / L(i, i);
end

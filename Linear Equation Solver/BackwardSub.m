%{
    BackwardSub(U, b) solves U * x = b, 
    where U is a upper triangular matrix. (Backward Substitution)
%}

function [x] = BackwardSub(U, b)

% Find size of L
row = length(U);

% Initialize solution
x = zeros(row, 1);

x(row) = b(row) / U(row, row);

% Backward substitution
for i = row - 1 : -1 : 1
    x(i) = b(i);
    for j = row : -1 : i + 1
        x(i) = x(i) - U(i, j) * x(j);
    end
    x(i) = x(i) / U(i, i);
end

function [x_next] = Jacobi(A, b, x, iter)

%{
    Jacobi Method solves the (square) system of linear equations A*x = b,
    with 'iter' as number of iterations
    x as initial guess
    Write A = D + R where, D := diag(a11, a22, ... ), R := A - D.
    The convergence condition is (spectral radius of D^-1 * R) < 1.
%}

% Find size of matrix A
row = length(A);

% Set initial guess of x as 1/n (1, 1, ... , 1)
x_prev = x;
x_next = zeros(row, 1);

% A = D + R, x_next = D^-1 ( b - R * x_prev)
for k = 1 : iter
    for i = 1 : row
        sum = 0;
        for j = 1 : row
            if i ~= j
                sum = sum + A(i, j) * x_prev(j, 1);
            end
        end
        x_next(i, 1) = (b(i, 1) - sum) / A(i, i);
    end
    x_prev = x_next;
end


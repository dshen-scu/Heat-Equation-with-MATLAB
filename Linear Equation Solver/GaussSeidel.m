function [x] = GaussSeidel(A, b, x, iter)

%{
    Gauss-Seidel Method,
    for solving a square system of linear equations A * x = b.
    with 'iter' as the number of iterations,
    x as initial guess for the solution.
    
    This method converges if A is:
    - symmetric positive definite
    - strictly or irreducibly diagonally dominant
%}

    row = length(A);
    for k = 1 : iter
        for i = 1 : row
            x(i) = (1 / A(i, i)) * (b(i) - A(i, 1 : row) * x + A(i, i) * x(i));
        end
    end
end
    
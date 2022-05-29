function x = solveEquationWithCholesky(A, b)
    A = cholesky(A);
    y = solveLowerTriangularEquation(A, b);
    x = solveUpperTriangularEquation(A.', y);
end
function x = solveEquationWithCholeskyImproved(A, b)
    [A, D] = choleskyImproved(A);
    y = solveLowerTriangularEquation(A, b)./(D.');
    x = solveUpperTriangularEquation(A.', y);
end
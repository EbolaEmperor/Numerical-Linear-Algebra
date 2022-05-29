function x = solveEquationWithLU(A, b)
    A = getLU(A);
    y = solveUnitLowerTriangularEquation(A, b);
    x = solveUpperTriangularEquation(A, y);
end


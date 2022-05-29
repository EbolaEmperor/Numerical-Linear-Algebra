function x = solveEquationWithPLU(A, b)
    [p,A] = getPLU(A);
    n = size(b,1);
    for i = 1:n-1
        b([i p(i)],:) = b([p(i) i],:);
    end
    y = solveUnitLowerTriangularEquation(A, b);
    x = solveUpperTriangularEquation(A, y);
end


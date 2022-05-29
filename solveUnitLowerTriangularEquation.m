function x = solveUnitLowerTriangularEquation(L, b)
    n = size(L,1);
    for i = 1:n-1
        b(i+1:n,1) = b(i+1:n,1) - b(i)*L(i+1:n,i);
    end
    x = b;
end
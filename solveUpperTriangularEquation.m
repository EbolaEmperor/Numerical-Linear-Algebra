function x = solveUpperTriangularEquation(U, b)
    n = size(U, 1);
    for i = n : -1 : 2
        b(i) = b(i) / U(i,i);
        b(1:i-1) = b(1:i-1) - b(i) * U(1:i-1, i);
    end
    b(1) = b(1) / U(1,1);
    x = b;
end
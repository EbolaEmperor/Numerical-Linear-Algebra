function x = solveLowerTriangularEquation(L, b)
    n = size(L,1);
    for i = 1:n-1
        b(i,1) = b(i,1)/L(i,i);
        b(i+1:n,1) = b(i+1:n,1) - b(i)*L(i+1:n,i);
    end
    b(n,1) = b(n,1)/L(n,n);
    x = b;
end
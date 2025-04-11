function A = genDiff(n)
    A = 2 * eye(n);
    A(1:n-1, 2:n) = A(1:n-1, 2:n) - eye(n-1);
    A(2:n, 1:n-1) = A(2:n, 1:n-1) - eye(n-1);
end
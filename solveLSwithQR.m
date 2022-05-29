function x = solveLSwithQR(A,b)
% 用QR分解求解最小二乘问题
    m = size(A,1);
    n = size(A,2);
    [Q,R] = getQR(A);
    c = Q(1:m,1:n).'*b;
    x = solveUpperTriangularEquation(R,c);
end


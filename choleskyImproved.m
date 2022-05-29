function [L, D] = choleskyImproved(A)
%   改进的cholesky分解算法，计算L和D使得A=LDL^T
    n = size(A,1);
    L = eye(n);
    D = zeros(1,n);
    for j = 1:n
        v = D(1:j-1).*L(j,1:j-1);
        D(j) = A(j,j) - dot(L(j,1:j-1),v);
        L(j+1:n,j) = (A(j+1:n,j)-L(j+1:n,1:j-1)*(v.'))/D(j);
    end
end
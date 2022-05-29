function L = cholesky(A)
%   cholesky分解算法，计算L使得A=LL^T
    n = size(A,1);
    L = zeros(n,n);
    for k = 1:n
        L(k,k) = sqrt(A(k,k) - sum(L(k,1:k-1).^2));
        L(k+1:n,k) = (A(k+1:n,k) - L(k+1:n,1:k-1)*(L(k,1:k-1).'))/L(k,k);
    end
end
function [A,Q] = tridiag(A)
% 计算矩阵的三对角分解
    n = size(A,1);
    Q = eye(n);
    for k = 1:n-2
        [v,beta] = householder(A(k+1:n,k));
        u = beta*(A(k+1:n,k+1:n)*v);
        w = u-(beta*(u'*v)/2)*v;
        A(k+1,k) = norm(A(k+1:n,k),2);
        A(k,k+1) = A(k+1,k);
        A(k+1:n,k+1:n) = A(k+1:n,k+1:n)-v*w'-w*v';
        A(k+2:n,k) = zeros(n-k-1,1);
        A(k,k+2:n) = zeros(1,n-k-1);
        Q = Q * blkdiag(eye(k),eye(n-k)-beta*(v*v'));
    end
end


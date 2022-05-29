function [H,Q] = hessenberg(A)
% 计算矩阵A的上Hessenberg分解，H=Q'AQ
    n = size(A,1);
    Q = eye(n);
    for k = 1:n-2
        [v,beta] = householder(A(k+1:n,k));
        Hk = (eye(n-k)-beta*(v*v'));
        A(k+1:n,k:n) = Hk*A(k+1:n,k:n);
        A(1:n,k+1:n) = A(1:n,k+1:n)*Hk;
        Q = Q * blkdiag(eye(k), Hk);
    end
    H = A;
    for k = 3:n
        H(k,1:k-2) = zeros(1,k-2);
    end
end


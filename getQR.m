function [Q,R] = getQR(A)
% 计算矩阵A的QR分解，其中A的行数不小于列数
    m = size(A,1);
    n = size(A,2);
    Q = eye(m);
    for j = 1:n
        if j<m
            [v,beta] = householder(A(j:m,j));
            H = eye(m-j+1)-beta*(v*v.');
            A(j:m,j:n) = H * A(j:m,j:n);
            Q = Q * blkdiag(eye(j-1), H);
        end
    end
    R = A(1:n,1:n);
end


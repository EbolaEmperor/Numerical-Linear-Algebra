function [T,Q] = symmetricEigen(A)
% 计算实对称矩阵的特征值和特征向量
    n = size(A,1);
    [T,Q] = tridiag(A);
    u = 1e-16;
    while true
        % 收敛性判定
        for i = 1:n-1
            if abs(T(i+1,i)) <= (abs(T(i,i))+abs(T(i+1,i+1)))*u
                T(i+1,i) = 0;
                T(i,i+1) = 0;
            end
        end
        m = 0;
        while m<n && (m==n-1 || abs(T(n-m,n-m-1))<u)
            m = m+1;
        end
        if m==n
            break;
        end
        l = n-m;
        while l>1 && abs(T(l,l-1))>=u
            l = l-1;
        end
        l = l-1;
        % 进入Wilkinson位移隐式对称QR迭代
        [T(l+1:n-m,l+1:n-m),G] = wilkinsonQR(T(l+1:n-m,l+1:n-m));
        Q = Q * blkdiag(eye(l),G',eye(m));
    end
    T = diag(diag(T));
end


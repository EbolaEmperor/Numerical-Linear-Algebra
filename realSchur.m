function [H,Q] = realSchur(A)
% REALSCHUR 计算实矩阵的实Schur分解：隐式QR算法
    n = size(A,1);
    [H,Q] = hessenberg(A);
    u = 1e-14;
    while true
        % 收敛性判定
        for i = 2:n
            if abs(H(i,i-1)) <= (abs(H(i,i))+abs(H(i-1,i-1)))*u
                H(i,i-1) = 0;
            end
        end
        m = 0;
        while m<n
            if m==n-1 || abs(H(n-m,n-m-1))<u
                m = m+1;
            else
                if (m==n-2 || abs(H(n-m-1,n-m-2))<u) && isComplexEigen(H(n-m-1:n-m,n-m-1:n-m))
                    m = m+2;
                else
                    break;
                end
            end
        end
        if m==n
            break;
        end
        l = n-m;
        while l>1 && abs(H(l,l-1))>=u
            l = l-1;
        end
        l = l-1;
        % 进入双重步隐式QR迭代
        [H(l+1:n-m,l+1:n-m),P] = doubleQR(H(l+1:n-m,l+1:n-m));
        Q = Q * blkdiag(eye(l),P,eye(m));
        H(1:l,l+1:n-m) = H(1:l,l+1:n-m)*P;
        H(l+1:n-m,n-m+1:n) = P'*H(l+1:n-m,n-m+1:n);
    end
    for k = 3:n
        H(k,1:k-2) = zeros(1,k-2);
    end
end


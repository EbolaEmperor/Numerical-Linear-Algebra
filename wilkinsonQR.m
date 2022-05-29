function [T,Q] = wilkinsonQR(T)
% 带Wilkinson位移的隐式对称QR迭代
    n = size(T,1);
    d = (T(n-1,n-1)-T(n,n))/2;
    mu = T(n,n)+d-sign(d)*sqrt(d*d+T(n,n-1)*T(n,n-1));
    x = T(1,1)-mu;
    z = T(2,1);
    Q = eye(n);
    for k = 1:n-1
        G = givensMatrix(x,z,k,k+1,n);
        T = G*T*G';
        Q = G*Q;
        if k < n-1
            x = T(k+1,k);
            z = T(k+2,k);
        end
    end
end


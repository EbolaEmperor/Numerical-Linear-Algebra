function [H,P] = doubleQR(H)
    % Francis双重步位移QR迭代算法
    n = size(H,1);
    if n <= 2
        mu = H(n, n);
        H = H - mu * eye(n);
        [P, R] = getQR(H);
        H = R * P + mu * eye(n);
        return;
    end
    P = eye(n);
    m = n-1;
    s = H(m,m)+H(n,n);
    t = H(m,m)*H(n,n)-H(m,n)*H(n,m);
    x = H(1,1)*H(1,1)+H(1,2)*H(2,1)-s*H(1,1)+t;
    y = H(2,1)*(H(1,1)+H(2,2)-s);
    if n >= 3
        z = H(2,1)*H(3,2);
    end
    for k = 0:n-3
        [v,beta] = householder([x,y,z]');
        q = max([1,k]);
        Pk = eye(3)-beta*(v*v');
        H(k+1:k+3,q:n) = Pk*H(k+1:k+3,q:n);
        r = min([k+4,n]);
        H(1:r,k+1:k+3) = H(1:r,k+1:k+3)*Pk;
        x = H(k+2,k+1);
        y = H(k+3,k+1);
        if k < n-3
            z = H(k+4,k+1);
        end
        P = P * blkdiag(eye(k),Pk,eye(n-k-3));
    end
    [v,beta] = householder([x,y]');
    Pk = eye(2)-beta*(v*v');
    if n >= 3
        H(n-1:n,n-2:n) = Pk*H(n-1:n,n-2:n);
    end
    H(1:n,n-1:n) = H(1:n,n-1:n)*Pk;
    P = P * blkdiag(eye(n-2),Pk);
    for k = 3:n
        H(k,1:k-2) = zeros(1,k-2);
    end
end


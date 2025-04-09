function [Q, R] = getQR(A)
    [m, n] = size(A);
    R = A;
    Q = eye(m);
    for k = 1:n
        x = R(k:m, k);
        [v, beta] = householder(x, 2.3e-16);
        if beta ~= 0
            H_k = eye(m-k+1) - beta * (v * v');
            R(k:m, k:n) = H_k * R(k:m, k:n);
            Q(k:m,:) = H_k * Q(k:m,:);
        end
    end
    R = triu(R);
    Q = Q';
end
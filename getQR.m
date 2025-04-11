function [Q, R] = getQR(A)
    [m, n] = size(A);
    R = A;
    Q = eye(m);
    for k = 1:n
        x = R(k:m, k);
        [v, beta] = householder(x, 2.3e-16);
        if beta ~= 0
            R(k:m, k:n) = R(k:m, k:n) - beta * v * (v' * R(k:m, k:n));
            Q(k:m,:) = Q(k:m,:) - beta * v * (v' * Q(k:m,:));
        end
    end
    R = triu(R);
    Q = Q';
end
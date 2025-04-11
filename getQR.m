function [Q, R] = getQR(R)
    [m, n] = size(R);
    Q = eye(m);
    for k = 1:n
        [v, beta] = householder(R(k:m, k));
        if beta ~= 0
            R(k:m, k:n) = R(k:m, k:n) - beta * v * (v' * R(k:m, k:n));
            Q(k:m,:) = Q(k:m,:) - beta * v * (v' * Q(k:m,:));
        end
    end
    R = triu(R);
    Q = Q';
end
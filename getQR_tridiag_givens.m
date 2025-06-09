function [Q, R] = getQR_tridiag_givens(R)
    [m, n] = size(R);
    Q = eye(m);
    for j = 1 : n-1
        i = j + 1;
        [c, s] = givens(R(i-1, j), R(i, j));
        
        r1 = R(i-1, j:n);
        r2 = R(i,   j:n);
        R(i-1, j:n) = c * r1 + s * r2;
        R(i,   j:n) = -s * r1 + c * r2;
        
        q1 = Q(i-1, :);
        q2 = Q(i,   :);
        Q(i-1, :) = c * q1 + s * q2;
        Q(i,   :) = -s * q1 + c * q2;
    end
    R = triu(R);
    Q = Q';
end
function [H, Q] = hessenberg(A)
% HESSENBERG 计算矩阵 A 的上 Hessenberg 分解，满足 H = Q' * A * Q
% 输入:
%   A - 待分解矩阵
% 输出:
%   H - 上 Hessenberg 矩阵（主对角线下仅有一条次对角线非零）
%   Q - 正交矩阵，使得 H = Q' * A * Q

    n = size(A, 1);
    Q = eye(n);
    H = A;
    for k = 1:n-2
        [v, beta] = householder(H(k+1:n, k));
        w = beta * (v' * H(k+1:n, k:n));
        H(k+1:n, k:n) = H(k+1:n, k:n) - v * w;
        w = beta * (H(1:n, k+1:n) * v);
        H(1:n, k+1:n) = H(1:n, k+1:n) - w * v';
        w = beta * (Q(:, k+1:n) * v);
        Q(:, k+1:n) = Q(:, k+1:n) - w * v';
    end
    H(tril(true(n), -2)) = 0;
end
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
        [v, beta] = householder(H(k+1:end, k));
        w = beta * (v' * H(k+1:end, k:end));
        H(k+1:end, k:end) = H(k+1:end, k:end) - v * w;
        w = beta * (H(1:end, k+1:end) * v);
        H(1:end, k+1:end) = H(1:end, k+1:end) - w * v';
        w = beta * (Q(:, k+1:end) * v);
        Q(:, k+1:end) = Q(:, k+1:end) - w * v';
    end
    H(tril(true(n), -2)) = 0;
end
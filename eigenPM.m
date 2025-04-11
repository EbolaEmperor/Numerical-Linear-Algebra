function [D, Q] = eigenPM(A, eps)
% 用幂法求矩阵的特征根和特征向量
    if nargin < 2
        eps = 2.3e-16;
    end
    n = size(A, 1);
    D = zeros(n, 1);
    Q = zeros(n, n);
    for i = 1 : n
        [D(i), Q(:,i)] = powerMethod(A, eps, 2000 * n);
        Q(:,i) = Q(:,i) / norm(Q(:,i));
        A = A - D(i) * Q(:,i) * Q(:,i)';
    end
end
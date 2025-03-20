function [x, step] = sor(A, b, omega, resTol, maxstep)
    if nargin < 5
        maxstep = 20000;
    end
    n = size(A, 1);
    x = b;
    step = 0;
    g = b ./ diag(A);
    B = eye(n);
    % 手动计算 D^{-1} * A，以避免 O(n^3) 的矩阵乘法
    for i = 1:n
        B(i,:) = B(i,:) - 1.0 / A(i,i) * A(i,:);
    end
    % 因为 matlab 里的稀疏矩阵采用列压缩存储，所以 B(:,i) 会比 B(i,:) 快得多
    % 而我们的 SOR 迭代需要频繁使用 B(i,:)，所以出于效率考虑，将 B 转置一下
    B = B';
    if nnz(A) < n^1.3
        A = sparse(A);
        B = sparse(B);
    end
    while norm(A * x - b) > resTol && step < maxstep
        for i = 1:n
            x(i) = x(i) * (1 - omega) + omega * g(i) + omega * dot(B(:,i), x);
        end
        step = step + 1;
        if mod(step, 100) == 0
            fprintf("Iteration %d\n", step);
        end
    end
end


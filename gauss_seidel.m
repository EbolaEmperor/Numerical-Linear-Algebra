function x = gauss_seidel(A, b, x0, tol, maxIter)
    n = length(b);
    x = x0;
    for iter = 1:maxIter
        x_old = x;
        for i = 1:n
            s1 = A(i,1:i-1) * x(1:i-1); % 左侧部分求和
            s2 = A(i,i+1:n) * x_old(i+1:n); % 右侧部分求和（使用旧解）
            x(i) = (b(i) - s1 - s2) / A(i,i); % 更新当前解
        end
        
        if max(abs(x - x_old)) <= tol
            fprintf('Converged in %d iterations.\n', iter);
            return;
        end
    end
    fprintf('Reach the maximum iteration times.\n');
end

function [D,Q,step] = symEigJacobi(A)
% 过关Jacobi方法求实对称矩阵的所有特征根和特征向量
    n = size(A,1);
    sigma = 2*n; %关值的递减系数
    delta = max(max(abs(A))) / 2;
    Q = eye(n);
    step = 0;
    while true
        pass = true;
        for p = 1:n
            for q = p+1:n
                if abs(A(p,q)) > delta
                    pass = false;
                    [c, s] = computeJacobiCoeffs(A, p, q);
                    
                    % 选出除 p,q 外的其它下标
                    I = setdiff(1:n, [p, q]);
                    
                    % 向量化更新 A(I, p) 和 A(I, q)
                    Ap_temp = A(I, p);
                    Aq_temp = A(I, q);
                    A(I, p) = c * Ap_temp - s * Aq_temp;
                    A(I, q) = s * Ap_temp + c * Aq_temp;
                    
                    % 由于 A 是对称矩阵，更新对应的行
                    A(p, I) = A(I, p)';
                    A(q, I) = A(I, q)';
                    
                    % 更新对角元素及 (p,q) 元
                    App = A(p, p);
                    Aqq = A(q, q);
                    Apq = A(p, q);
                    A(p, p) = c^2 * App - 2 * c * s * Apq + s^2 * Aqq;
                    A(q, q) = s^2 * App + 2 * c * s * Apq + c^2 * Aqq;
                    A(p, q) = 0;
                    A(q, p) = 0;
                    
                    % 向量化更新特征向量矩阵 Q
                    temp = Q(:, p);
                    Q(:, p) = c * Q(:, p) - s * Q(:, q);
                    Q(:, q) = s * temp + c * Q(:, q);
                end
            end
        end
        if pass, break; end
        step = step+1;
        delta = delta / sigma;
    end
    D = diag(diag(A));
end

function [c, s] = computeJacobiCoeffs(A, p, q)
    % 计算 Jacobi 旋转的系数，对于下标 p 和 q
    %
    % 公式:
    % \(
    % \tau = \frac{A(q,q)-A(p,p)}{2A(p,q)}
    % \)
    % 若 \(\tau \ge 0\)，则取
    % \(
    % t = \frac{1}{\tau+\sqrt{1+\tau^2}}
    % \)
    % 否则取
    % \(
    % t = -\frac{1}{-\tau+\sqrt{1+\tau^2}}
    % \)
    % 然后有
    % \(
    % c = \frac{1}{\sqrt{1+t^2}},\quad s = t\,c.
    % \)
    if A(p, q) == 0
        c = 1;
        s = 0;
    else
        tau = (A(q, q) - A(p, p)) / (2 * A(p, q));
        if tau >= 0
            t = 1 / (tau + sqrt(1 + tau^2));
        else
            t = -1 / (-tau + sqrt(1 + tau^2));
        end
        c = 1 / sqrt(1 + t^2);
        s = t * c;
    end
end
function [D,Q,step] = symEigJacobi(A)
% 过关Jacobi方法求实对称矩阵的所有特征根和特征向量
    eps = 1e-16;
    n = size(A,1);
    sigma = 2*n; %关值的递减系数
    delta = max(max(abs(A)));
    Q = eye(n);
    step = 0;
    while true
        pass = true;
        for p = 1:n
            for q = p+1:n
                if abs(A(p,q))>eps
                    P = symShar2(A,p,q);
                    A = P'*A*P;
                    Q = Q*P;
                    pass = false;
                end
            end
        end
        if pass, break; end
        step = step+1;
        delta = delta / sigma;
    end
    D = diag(diag(A));
end

function Q = symShar2(A,p,q)
    tau = (A(q,q)-A(p,p))/(2*A(p,q));
    if tau==0, t=1;
    else, t=sign(tau)/(abs(tau)+sqrt(1+tau^2)); end
    c = 1/sqrt(1+t^2);
    s = t*c;
    Q = eye(size(A,1));
    Q(p,p) = c; Q(p,q) = s;
    Q(q,p) = -s; Q(q,q) = c;
end
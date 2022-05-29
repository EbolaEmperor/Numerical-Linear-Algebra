function [v,beta] = householder(x)
% 计算x的Householder变换
    n = length(x);
    v = zeros(n,1);
    if norm(x,inf)==0
        beta = 0;
        return;
    end
    ita = norm(x,inf);
    x = x/ita;
    sigma = x(2:n).'*x(2:n);
    v(2:n) = x(2:n);
    if sigma==0
        beta = 0;
    else
        alpha = sqrt(x(1)*x(1)+sigma);
        if x(1)<=0
            v(1) = x(1)-alpha;
        else
            v(1) = -sigma/(x(1)+alpha);
        end
        beta = 2*v(1)*v(1)/(sigma+v(1)*v(1));
        v = v/v(1);
    end
end


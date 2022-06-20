function [lambda,v] = dichotomy(T,m)
% 利用二分法求非退化对称三对角阵T的第m大特征值
    r = norm(T,inf);
    l = -r;
    eps = 1e-15;
    while r-l > eps
        mid = (l+r)/2;
        if variant(T,mid) >= m, r=mid;
        else, l=mid; end
    end
    lambda = l;
    n = size(T,1);
    v = inversePM(T,lambda,ones(n,1));
end

function m = variant(T,x)
% 求变号数
    m = 0;
    q = T(1,1) - x;
    if q<0, m=m+1; end
    n = size(T,1);
    for i = 2:n
        q = T(i,i) - x - T(i,i-1)^2/q;
        if q<0, m=m+1; end
    end
end
function [root] = maxroot(a)
    n = size(a,1);
    A = zeros(n,n);
    A(2:n,1:n-1) = eye(n-1);
    A(1:n,n) = -a;
    root = maxeig(A);
end


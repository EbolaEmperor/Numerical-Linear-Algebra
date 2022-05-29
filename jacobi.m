function [x,step] = jacobi(A,b,err)
    n = size(A,1);
    H = diag(ones(n,1)./diag(A));
    B = eye(n) - H*A;
    g = H*b;
    x0 = zeros(n,1);
    x = g;
    step = 0;
    while vecnorm(x-x0)>err
        x0 = x;
        x = B*x + g;
        step = step + 1;
    end
end


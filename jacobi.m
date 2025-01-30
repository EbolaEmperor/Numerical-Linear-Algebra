function [x,step] = jacobi(A,b,err,x0)
    n = size(A,1);
    H = diag(ones(n,1)./diag(A));
    B = eye(n) - H*A;
    g = H*b;
    if nargin < 4
        x = zeros(n,1);
    else
        x = x0;
    end
    step = 0;
    while step == 0 || vecnorm(x-x0)>err
        x0 = x;
        x = B*x + g;
        step = step + 1;
    end
end


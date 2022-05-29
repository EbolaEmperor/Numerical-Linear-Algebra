function [x,step] = sor(A,b,omega,err)
    n = size(A,1);
    x0 = zeros(n,1);
    x = b;
    step = 0;
    g = b./diag(A);
    B = eye(n) - diag(ones(n,1)./diag(A))*A;
    while vecnorm(x-x0)>err && step<20000
        x0 = x;
        for i = 1:n
            x(i) = x(i) * (1-omega) + omega*g(i);
            for j = 1:n
                if j~=i
                    x(i) = x(i) + omega*B(i,j)*x(j);
                end
            end
        end
        step = step + 1;
    end
end


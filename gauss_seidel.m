function [x,step] = gauss_seidel(A,b,err)
    n = size(A,1);
    x0 = zeros(n,1);
    x = b;
    step = 0;
    while vecnorm(x-x0)>err
        x0 = x;
        for i = 1:n
            x(i) = b(i);
            for j = 1:n
                if j~=i
                    x(i) = x(i) - A(i,j)*x(j);
                end
            end
            x(i) = x(i)/A(i,i);
        end
        step = step + 1;
    end
end


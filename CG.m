function [x,step] = CG(A,b,err)
    x = zeros(size(A,1),1);
    r = A*x-b;
    p = -r;
    step = 0;
    while vecnorm(r) >= err
        step = step + 1;
        alpha = - (r.'*p)/(p.'*A*p);
        x = x + alpha*p;
        tmp = r.'*r;
        r = r + alpha*A*p;
        beta = (r.'*r)/tmp;
        p = -r + beta*p;
    end
end


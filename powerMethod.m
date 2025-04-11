function [lambda, q] = powerMethod(A, eps, times, q)
    if nargin < 4
        q = rand(size(A,2), 1);
    end
    if nargin < 3
        times = size(A,1) * 1000;
    end
    if nargin < 2
        eps = 2.3e-16;
    end
    lambda = 1;
    for k = 1:times
        z = A * q;
        [~, r] = max(abs(z));
        if(abs((z(r)-lambda)/z(r)) <= eps)
            fprintf("Terminated at iter %d\n", k);
            return;
        end
        q = z / z(r);
        lambda = z(r);
    end
    fprintf("[Warning] Power method did not converge after iter %d\n", times);
end
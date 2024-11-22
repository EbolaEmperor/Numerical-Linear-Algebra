function [lambda, q] = powerMethod(A, eps, times)
    q = rand(size(A,2), 1);
    lambda = 1;
    for k = 1:times
        z = A * q;
        if(norm(z - lambda * q) < eps)
            fprintf("Terminated at iter %d\n", k);
            return;
        end
        q = z / norm(z);
        lambda = q' * A * q;
    end
    fprintf("[Warning] Power method did not converge after iter %d\n", times);
end
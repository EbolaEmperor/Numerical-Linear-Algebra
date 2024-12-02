function [lambda, q] = RayleighIter(A, eps, times)
    n = size(A, 1);
    q = rand(n, 1);
    q = q / norm(q);
    lambda = 1;
    for k = 1:times
        lambda = q' * A * q;
        q = (A - lambda * eye(n)) \ q;
        q = q / norm(q);
        if(norm(A * q - lambda * q) < eps)
            fprintf("Terminated at iter %d\n", k);
            return;
        end
    end
    fprintf("[Warning] Rayleigh quotient iteration did not converge after iter %d\n", times);
end
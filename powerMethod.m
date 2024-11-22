function [lambda, q] = powerMethod(A, eps, times)
    q = rand(size(A,2), 1);
    lambda = 1;
    oldlam = 1;
    cnt = 0;
    for k = 1:times
        z = A * q;
        if(norm(z - lambda * q) < eps)
            fprintf("Terminated at iter %d\n", k);
            return;
        end
        q = z / norm(z);
        oldoldlam = oldlam;
        oldlam = lambda;
        lambda = q' * A * q;
        rate = abs((oldoldlam-oldlam) / (oldlam-lambda));
        if k >= 3
            fprintf("Iter %d: reduction rate %f\n", k, rate);
            if rate < 1.1
                cnt = cnt + 1;
            else
                cnt = 0;
            end
            if cnt > 5
                fprintf("[Warning] Terminated at iter %d because of no more reduction\n", k);
                return;
            end
        end
    end
    fprintf("[Warning] Power method did not converge after iter %d\n", times);
end
function [v, beta] = householder(x, eps)
    if nargin < 2
        eps = 0;
    end
    if isscalar(x)
        beta = 0;
        v = x;
        return;
    end
    sigma = dot(x(2:end), x(2:end));
    v = x;
    if sigma <= eps
        beta = 0;
    else
        mu = sqrt(x(1)^2 + sigma);
        if x(1) <= 0
            v(1) = x(1) - mu;
        else
            v(1) = -sigma / (x(1) + mu);
        end
        beta = 2 * (v(1)^2) / (sigma + v(1)^2);
        v = v / v(1);
    end
end
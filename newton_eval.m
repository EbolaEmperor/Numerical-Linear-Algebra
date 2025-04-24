function p = newton_eval(x, a, z)
    n = numel(a);
    z = z(:).';
    p = a(n) * ones(size(z));
    for k = n-1:-1:1
        p = a(k) + (z - x(k)) .* p;
    end
end
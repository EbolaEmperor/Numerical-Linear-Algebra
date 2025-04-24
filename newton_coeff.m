function a = newton_coeff(x, y)
    n = numel(x);
    a = y(:);
    for j = 2:n
        a(j:n) = (a(j:n) - a(j-1:n-1)) ./ (x(j:n) - x(1:n-j+1)).';
    end
end
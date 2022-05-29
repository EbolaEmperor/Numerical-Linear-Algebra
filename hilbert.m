function [A] = hilbert(n)
    A = zeros(n,n);
    for i = 1:n
        for j = 1:n
            A(i,j) = 1.0/(i+j-1);
        end
    end
end

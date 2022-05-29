function [p, A] = getPLU(A)
%   列主元法高斯消元计算LU分解
    n = size(A,1);
    p = zeros(1,n);
    for j = 1:n-1
        p(j) = j;
        for k = j+1:n
            if abs(A(k,j))>abs(A(p(j),j))
                p(j) = k;
            end
        end
        A([j p(j)],:) = A([p(j),j],:);
        A(j+1:n,j) = A(j+1:n,j)/A(j,j);
        A(j+1:n,j+1:n) = A(j+1:n,j+1:n) - A(j+1:n,j)*A(j,j+1:n);
    end
end


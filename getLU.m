function A = getLU(A)
%   高斯消元法计算LU分解
    n = size(A,1);
    for j=1:n-1
        A(j+1:n,j) = A(j+1:n,j)/A(j,j);
        A(j+1:n,j+1:n) = A(j+1:n,j+1:n) - A(j+1:n,j)*A(j,j+1:n);
    end
end
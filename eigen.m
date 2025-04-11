function D = eigen(A)
% 求矩阵的特征根
    n = size(A,1);
    [H,~] = realSchur(A);
    D = zeros(n,1);
    i = 1;
    while i<=n
        if i==n || H(i+1,i)==0
            D(i) = H(i,i);
            i = i+1;
        else
            D(i:i+1) = getComplexEigen(H(i:i+1,i:i+1));
            i = i+2;
        end
    end
    D = sort(D, "descend");
end

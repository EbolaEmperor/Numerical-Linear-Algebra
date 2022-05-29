function G = givensMatrix(a,b,i,j,n)
% 返回对i,j行列变换的givens变换矩阵
    [c,s] = givens(a,b);
    G = eye(n);
    G(i,i) = c;
    G(i,j) = s;
    G(j,i) = -s;
    G(j,j) = c;
end


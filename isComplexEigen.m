function res = isComplexEigen(A)
% 判断一个2*2实矩阵的特征值是不是复数
    res = ((A(1,1)+A(2,2))*(A(1,1)+A(2,2)) - 4*det(A) < 0);
end


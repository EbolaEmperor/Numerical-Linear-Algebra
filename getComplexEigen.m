function res = getComplexEigen(A)
% 求一个2*2实矩阵的复特征值
    res = zeros(2,1);
    delta = (A(1,1)+A(2,2))*(A(1,1)+A(2,2)) - 4*det(A);
    res(1) = (A(1,1)+A(2,2)+sqrt(delta))/2;
    res(2) = (A(1,1)+A(2,2)-sqrt(delta))/2;
end


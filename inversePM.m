function v = inversePM(A,lambda,z)
% 反幂法计算特征向量
    n = size(A,1);
    y = (A-lambda*eye(n))\z;
    v = y/vecnorm(y);
end


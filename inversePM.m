function v = inversePM(A,lambda,z)
% 反幂法计算特征向量
    n = size(A,1);
    eps = 1e-15;
    step = 1;
    v = z;
    while vecnorm(A*v-lambda*v)>=eps && step <= 2
        y = (A-lambda*eye(n))\v;
        v = y/vecnorm(y);
        step = step + 1;
    end
end


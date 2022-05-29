function lambda = maxeig(A)
    x = rand(size(A,1),1);
    step = 0;
    while(step<5000)
        step = step + 1;
        x1 = A*x;
        y1 = x1 / norm(x1,inf);
        lambda = max(abs(x1));
        % 若序列收敛，则有正的实特征根
        if(norm(y1-x,inf)<1e-5)
            return;
        end
        x = y1;
    end
    x = x / norm(x,inf);
    x1 = A*x;
    y1 = x1 / norm(x1,inf);
    
    % 判断是否为负的实特征根
    if(vecnorm(x+y1)<1e-4)
        lambda = -lambda;
        return;
    end
    
    % 判断是否为两个互为相反数的实特征根
    x2 = A*x1;
    lambda1 = sqrt(norm(x2,inf));
    lambda = [lambda1; -lambda1];
    x = [x2+lambda1*x1, x2-lambda1*x1];
    if(vecnorm(A*x(:,1)-lambda1*x(:,1))<1e-4 && vecnorm(A*x(:,2)+lambda1*x(:,2))<1e-4)
        return;
    end
    
    % 判断是否为两个共轭复特征根
    x3 = A*x2;
    B = [x2(1:2), x1(1:2)];
    b = -x3(1:2);
    sol = B\b;
    p = sol(1);
    q = sol(2);
    if(vecnorm(x3+p*x2+q*x1)<1e-4)
        lambda = [complex(p/2,sqrt(q-p*p/4)); complex(p/2,-sqrt(q-p*p/4))];
        return;
    end
    
    lambda = [];
end


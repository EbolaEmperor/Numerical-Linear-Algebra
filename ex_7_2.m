n = 100;
A = diag(repmat([2], 1, n))+diag(repmat([-1], 1, n-1), 1)+diag(repmat([-1], 1, n-1), -1);
[lambda_min, v_min] = dichotomy(A,1)
err_min = vecnorm(A*v_min-lambda_min*v_min)
[lambda_max, v_max] = dichotomy(A,n)
err_max = vecnorm(A*v_max-lambda_max*v_max)

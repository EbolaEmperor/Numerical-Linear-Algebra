function A = iterQR(A, eps, iterNum)
    if nargin < 3, iterNum = 100 * size(A, 1); end
    if nargin < 2, eps = 1e-6; end
    tic
    A = hessenberg(A);
    time1 = toc;
    for iter = 1:iterNum
        [Q, R] = getQR(A);
        f = diag(A);
        A = R * Q;
        err = norm(diag(A) - f, inf);
        fprintf("Iteration %d, err = %e\n", iter, err);
        if err < eps, break; end
    end
    time2 = toc;
    fprintf("CPU time: %f sec (Hessenberg: %f sec, QR-iteration: %f sec)\n", time2, time1, time2-time1);
end
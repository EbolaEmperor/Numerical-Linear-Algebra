function x = solveLaplacian(n, b)
    A = genDiff2(n);
    % omega = 1.98;
    % D = spdiags(diag(A), 0, size(A,1), size(A,1));
    % M1 = 4 * (D/omega + tril(A, -1));
    % M2 = (omega / (2 - omega)) * (D/omega + triu(A, 1));
    solver = Lap2DMGSolver;
    Mfun = @(b) solver.FMGCycle(b, n);
    
    tol = 5e-6;
    maxit = 1000;
    x = pcg(A, b, tol, maxit, Mfun);
end

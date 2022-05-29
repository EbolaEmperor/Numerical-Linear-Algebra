function a = kappa(A,k)
    if k==1
        a = norm1Optimize(A)*norm1Optimize(inv(A));
    end
    if k==inf
        a = norm1Optimize(A.')*norm1Optimize(inv(A.'));
    end
end
function [L,d,E] = cholesky_GM(G)
%CHOLESKY_GM Gill-Murray 修正Cholesky分解
    n = size(G,1);
    gamma = max(abs(diag(G)));
    xi = max(abs(G-diag(diag(G))));
    nu = max([1, sqrt(n*n-1)]);
    beta2 = max([gamma, xi/nu, 1e-6]);
    c = diag(diag(G));
    L = eye(n);
    d = zeros(1,n);
    E = zeros(1,n);
    for j = 1:n
        q = j;
        for k = j+1:n
            if abs(c(k,k))>abs(c(q,q))
                q = k;
            end
        end
        G([j q],:) = G([q j],:);
        G(:,[j q]) = G(:,[q j]);
        L(j,1:j-1) = c(j,1:j-1)./d(1:j-1);
        c(j+1:n,j) = G(j+1:n,j) - c(j+1:n,1:j-1)*(L(j,1:j-1).');
        if j==n
            theta = 0;
        else
            theta = max(abs(c(j+1:n,j)));
        end
        d(j) = max([1e-3, abs(c(j,j)), theta*theta/beta2]);
        E(j) = d(j) - c(j,j);
        for i = j+1:n
            c(i,i) = c(i,i) - c(i,j)*c(i,j)/d(j);
        end
    end
end


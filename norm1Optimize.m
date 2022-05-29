function [nrm] = norm1Optimize(B)
    k = 1;
    n = size(B,1);
    x = (1.0/n)*ones(n,1);
    while k==1
        w = B*x;
        v = sign(w);
        z = (B.')*v;
        nmz = norm(z,inf);
        if nmz <= (z.')*x
            nrm = norm(w,1);
            k = 0;
        else
            x = zeros(n,1);
            for j = 1:n
                if abs(z(j))==nmz
                    x(j) = 1;
                    break;
                end
            end
            k = 1;
        end
    end
end

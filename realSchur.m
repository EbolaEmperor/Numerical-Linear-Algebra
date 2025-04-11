function [H, Q] = realSchur(A)
% REALSCHUR 计算实矩阵的实 Schur 分解：隐式 QR 算法
% 输入:
%   A - 实矩阵
% 输出:
%   H - 实 Schur 形式（上 quasi-triangular）
%   Q - 正交矩阵，使得 A = Q*H*Q'

    n = size(A,1);
    [H, Q] = hessenberg(A);
    u = 1e-14;
    while true
        d = abs(diag(H));
        sub = abs(diag(H,-1));
        threshold = (d(1:end-1) + d(2:end)) * u;
        idx = find(sub <= threshold);
        if ~isempty(idx)
            H(sub2ind(size(H), idx+1, idx)) = 0;
        end
        m = 0;
        while m < n
            if m == n-1 || abs(H(n-m, n-m-1)) < u
                m = m + 1;
            else
                if (m == n-2 || abs(H(n-m-1, n-m-2)) < u) && isComplexEigen(H(n-m-1:n-m, n-m-1:n-m))
                    m = m + 2;
                else
                    break;
                end
            end
        end
        if m == n
            break;
        end
        l = n - m;
        while l > 1 && abs(H(l, l-1)) >= u
            l = l - 1;
        end
        l = l - 1;
        %% 进入双重步隐式 QR 迭代
        block = l+1 : n-m;
        [H(block, block), P] = doubleQR(H(block, block));
        if ~isempty(block)
            Q(:, block) = Q(:, block) * P;
        end
        if l >= 1
            H(1:l, block) = H(1:l, block) * P;
        end
        if n-m < n
            H(block, n-m+1:n) = P' * H(block, n-m+1:n);
        end
    end
    H(tril(true(n), -2)) = 0;
end
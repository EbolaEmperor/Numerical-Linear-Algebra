classdef Lap2DMGSolver
    methods (Access = private)
        %% 构造二维拉普拉斯矩阵（向量化版）
        % 对于网格大小 \(n\)，构造矩阵 \(A\in\mathbb{R}^{n^2\times n^2}\) 满足
        % 主对角元为 \(4\)，相邻点对应系数为 \(-1\)。
        function A = getLaplacian2D(obj, n)
            % 利用 Kronecker 积构造二维离散拉普拉斯算子
            e = ones(n,1);
            T = spdiags([-e, 4*e, -e], -1:1, n, n);
            S = spdiags([-e, -e], [-1,1], n, n);
            A = kron(speye(n), T) + kron(S, speye(n));
            A = full(A);  % 若需要稠密矩阵则转换为 full
        end
        
        %% 限制算子（向量化）
        % 将细网格向量 \(x\)（对应 \(n\times n\) 网格）限制到粗网格，
        % 利用每个粗网格点对应的 \(2\times2\) 块的和：
        % \(
        % y_{ij} = x_{2i-1,2j-1}+ x_{2i-1,2j}+ x_{2i,2j-1}+ x_{2i,2j}
        % \)
        function y = restriction(obj, x, n)
            U = reshape(x, n, n);
            Y = U(1:2:end, 1:2:end) + U(1:2:end, 2:2:end) + ...
                U(2:2:end, 1:2:end) + U(2:2:end, 2:2:end);
            y = Y(:);
        end
        
        %% 延拓算子（向量化）
        % 将粗网格向量 \(x\)（对应 \(n\times n\) 网格）延拓到细网格，
        % 细网格尺寸为 \(2n\times 2n\)，直接复制对应点的值，
        % 即利用函数 \(\texttt{repelem}\) 完成操作：
        % \(
        % y_{ij} = x_{\lfloor (i+1)/2\rfloor,\lfloor (j+1)/2\rfloor}
        % \)
        function y = prolongation(obj, x, n)
            Y = reshape(x, n, n);
            U = repelem(Y, 2, 2);
            y = U(:);
        end
        
        %% 线性延拓（向量化实现）
        % 采用加权平均将粗网格向量 \(x\) 延拓到细网格，
        % 权重分布为：
        % 中心：\(4/9\)，邻边：\(2/9\)，角落：\(1/9\)。
        % 这里先将粗网格数据插入到细网格的奇数位置，然后利用卷积
        % 与核
        % \(
        % K=\begin{pmatrix} 1/9 & 2/9 & 1/9\\ 2/9 & 4/9 & 2/9\\ 1/9 & 2/9 & 1/9 \end{pmatrix}
        % \)
        % 进行扩展。
        function y = prolongationLinear(obj, x, n)
            m = 2 * n;
            Y = reshape(x, n, n);
            U = zeros(m, m);
            U(1:2:end, 1:2:end) = Y;
            K = [1/9, 2/9, 1/9; 2/9, 4/9, 2/9; 1/9, 2/9, 1/9];
            V = conv2(U, K, 'same');
            y = V(:);
        end
        
        %% 加权雅可比迭代（向量化）
        % 更新公式为
        % \(
        % x_{ij}^{\text{new}}=(1-w)x_{ij}+\frac{w}{4}\Big(b_{ij}+\sum_{\text{邻域}}x\Big)
        % \)
        function new_x = wJacobi(obj, x, b, w)
            if nargin < 4, w = 0.8; end
            n = sqrt(length(x));
            U = reshape(x, n, n);
            B = reshape(b, n, n);
            % 累加上下左右邻域（边界自动忽略不存在的部分）
            sum_val = zeros(n, n);
            sum_val(2:end, :)   = sum_val(2:end, :)   + U(1:end-1, :);
            sum_val(1:end-1, :) = sum_val(1:end-1, :) + U(2:end, :);
            sum_val(:, 2:end)   = sum_val(:, 2:end)   + U(:, 1:end-1);
            sum_val(:, 1:end-1) = sum_val(:, 1:end-1) + U(:, 2:end);
            newU = (1-w)*U + w*(B + sum_val)/4;
            new_x = newU(:);
        end
        
        %% 矩阵向量乘积（向量化）
        % 计算二维拉普拉斯算子 \(A\) 作用于向量 \(x\) 后的结果，
        % 即 \(\, y = Ax \, \)。
        function y = vmult(obj, x)
            n = sqrt(length(x));
            U = reshape(x, n, n);
            % 采用五点差分格式实现：  
            % \(
            % y_{ij} = 4U_{ij} - U_{i-1,j} - U_{i+1,j} - U_{i,j-1} - U_{i,j+1}
            % \)
            ymat = 4*U;
            ymat(2:end, :)   = ymat(2:end, :)   - U(1:end-1, :);
            ymat(1:end-1, :) = ymat(1:end-1, :) - U(2:end, :);
            ymat(:, 2:end)   = ymat(:, 2:end)   - U(:, 1:end-1);
            ymat(:, 1:end-1) = ymat(:, 1:end-1) - U(:, 2:end);
            y = ymat(:);
        end
    end
    
    methods (Access = public)
        %% V-周期多重网格方法
        % 先在细网格上进行 4 次加权雅可比平滑，再计算残差
        % \(\, r = b - Ax\,\)，限制到粗网格上递归求解，再延拓回细网格并平滑。
        function x = VCycle(obj, x, b, n)
            if n <= 4
                A = obj.getLaplacian2D(n);
                x = A \ b;
                return;
            end
            for k = 1:4
                x = obj.wJacobi(x, b);
            end
            r = b - obj.vmult(x);
            r_coarse = obj.restriction(r, n);
            e2 = obj.VCycle(zeros((n/2)*(n/2),1), r_coarse, n/2);
            x = x + obj.prolongation(e2, n/2);
            for k = 1:4
                x = obj.wJacobi(x, b);
            end
        end
        
        %% FMG 循环
        % 先在粗网格上计算近似解，再延拓到细网格上进行 VCycle 修正。
        function x = FMGCycle(obj, b, n)
            if n <= 4
                A = obj.getLaplacian2D(n);
                x = A \ b;
                return;
            end
            x0 = obj.FMGCycle(obj.restriction(b, n), n/2);
            x = obj.VCycle(obj.prolongation(x0, n/2), b, n);
        end
    end
end
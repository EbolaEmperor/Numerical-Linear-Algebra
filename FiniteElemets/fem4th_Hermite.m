function [u, d, K] = fem4th_Hermite(n, f)
    % fem4th_Hermite 利用 3 阶 Hermite 元求解四阶边值问题
    %   求解问题: u''''(x) = f(x),  x in (0,1),
    %            u(0)=u(1)=u'(0)=u'(1)=0.
    %
    % 输入参数:
    %   n : 将 [0,1] 均匀剖分为 n 个单元（共 n+1 个节点）
    %
    % 示例调用:
    %   >> [u, d] = fem4th_Hermite(10, @(x) 1.0);
    %   u 返回数值解，d 返回数值解的梯度
    
    if nargin < 1
        n = 10;  % 默认剖分成 10 个单元
    end
    
    % 网格剖分
    h = 1/n;
    x_nodes = linspace(0,1,n+1)';  % 节点坐标
    
    % 总自由度: 每个节点有两个自由度 (位移和斜率)
    ndof = 2*(n+1);
    
    % 初始化全局刚度矩阵 K 和载荷向量 F
    K = zeros(ndof, ndof);
    F = zeros(ndof, 1);
    
    % 3 点 Gauss-Legendre 求积 (在参考区间 [0,1] 上)
    % 参考区间 [0,1] 的求积点与权重:
    qp = [0.5*(1 - sqrt(3/5)), 0.5, 0.5*(1 + sqrt(3/5))];
    wq = [5/18, 4/9, 5/18];
    
    % 单元循环：每个单元连接节点 i 和 i+1
    for e = 1:n
        % 当前单元左端节点编号
        i = e;
        % 全局自由度对应顺序: [u_i, u'_i, u_{i+1}, u'_{i+1}]
        gdof = [2*i-1, 2*i, 2*(i+1)-1, 2*(i+1)];
        
        he = h;  % 单元长度（均匀剖分）
        
        % 单元刚度矩阵（参考上式）
        Ke = (1/he^3)*[12,      6*he,   -12,      6*he;
                       6*he,   4*he^2,  -6*he,    2*he^2;
                       -12,    -6*he,    12,     -6*he;
                       6*he,   2*he^2,  -6*he,    4*he^2];
                   
        % 组装到全局刚度矩阵中
        K(gdof, gdof) = K(gdof, gdof) + Ke;
        
        % 单元载荷向量初始化
        Fe = zeros(4,1);
        
        % 在参考区间 [0,1] 上积分，转换关系: x = x_nodes(i) + he*xi, dx = he*dxi
        for q = 1:length(qp)
            xi = qp(q);
            w  = wq(q);
            
            % 物理坐标
            xq = x_nodes(i) + he*xi;
            f_val = f(xq);
            
            % Hermite 形函数定义 (注意：对于 u' 相关形函数已乘以 he)
            N1 = 1 - 3*xi^2 + 2*xi^3;
            N2 = he*(xi - 2*xi^2 + xi^3);
            N3 = 3*xi^2 - 2*xi^3;
            N4 = he*(-xi^2 + xi^3);
            N  = [N1; N2; N3; N4];
            
            % 累加单元载荷 (积分时记得乘以雅可比因子 he)
            Fe = Fe + f_val * N * he * w;
        end
        
        % 组装载荷向量
        F(gdof) = F(gdof) + Fe;
    end
    
    % 施加边界条件: u(0)=0, u'(0)=0, u(1)=0, u'(1)=0
    % 对应全局自由度: 1, 2, ndof-1, ndof
    fixed_dofs = [1, 2, ndof-1, ndof];
    for k = fixed_dofs
        K(k,:) = 0;
        K(:,k) = 0;
        K(k,k) = 1;
        F(k) = 0;
    end
    
    % 求解线性方程组
    d = K \ F;
    
    % 提取节点处位移 (自由度序号 1,3,5,... 为 u)
    u = d(1:2:end);
    d = d(2:2:end);
    
    % % 在图中画出数值解
    % set(gcf,'Units','centimeters','Position',[6 6 50 15]);
    % subplot(1, 2, 1)
    % plot(x_nodes, u, '-o', 'LineWidth',1.5);
    % title('$u(x)$', 'Interpreter', 'latex');
    % grid on;
    % 
    % subplot(1, 2, 2)
    % plot(x_nodes, d, '-o', 'LineWidth',1.5);
    % title("$u'(x)$", 'Interpreter', 'latex');
    % grid on;
end
function [P_cushion_Pa, P_internal_state] = calc_cushion_pressure(N_FAN, S_leak_ft2, dzdt_ft_s, P_guess_psf)
% CALC_CUSHION_PRESSURE 计算气垫船气室压强 (4风机 -> 4气室版)
%
% 输入参数:
%   N_FAN       : [1x4] 向量, 四个风机的转速 [FR, FL, RR, RL] (RPM)
%   S_leak_ft2  : [1x4] 向量, 4个气室的底部泄流面积 (平方英尺)
%   dzdt_ft_s   : [1x4] 向量, 船体各气室处的垂直升沉速度 (英尺/秒) 
%   P_guess_psf : [4x1] 向量, 压强初值猜测 (psf) - 注意这里变成了4维
%
% 输出参数:
%   P_cushion_Pa     : [1x4] 向量, 4个气室的最终压强 (帕斯卡 Pa)
%   P_internal_state : [4x1] 向量,用于下一次迭代的初值 (psf)

    %% 常数定义 
    Area_cushion_ft2 = 800; % 单个气室面积
    psf_to_pa = 47.8803;    % 单位转换系数
    
    % 求解器配置
    max_iter = 30;
    tol = 20; % 容差
    lambda = 0.8; % 牛顿法步长因子 (简化模型后可以稍微调大一点)

    %% 检查输入
    % 注意：状态量现在只有4个 (4个气室)，不需要管道压强了
    if nargin < 4 || isempty(P_guess_psf) || length(P_guess_psf) ~= 4
        P = ones(4,1) * 109; % 给一个接近平衡状态的初值(例如109 psf)
    else
        P = P_guess_psf;
    end
    
    % 确保输入是列向量
    N_vec = N_FAN(:); 
    S = S_leak_ft2(:); 
    dzdt_vec = dzdt_ft_s(:);

    if length(N_vec) ~= 4
        error('错误：calc_cushion_pressure 需要 4 个风机转速输入');
    end

    %% 牛顿迭代求解
    converged = false;
    
    for k = 1:max_iter
        % 计算残差 F (现在是 4x1 向量)
        F = residual_function(P, S, N_vec, dzdt_vec, Area_cushion_ft2);
        
        % 检查收敛
        if norm(F) < tol
            converged = true;
            break; 
        end
        
        % 计算雅可比矩阵 (现在是 4x4 矩阵)
        J = compute_jacobian(P, S, N_vec, dzdt_vec, Area_cushion_ft2);
        
        % 更新步长
        delta = -J \ F;
        P = P + lambda * delta;
        
        P = max(P, 0.1); % 防止压强为负
    end
    
    if ~converged
        % 仿真初期容易不收敛，用警告代替报错，防止程序中断
        % warning('Pressure solver did not converge.');
    end

    %% 输出处理
    P_internal_state = P; % 保持 psf 单位
    
    % 直接输出这4个压强
    P_cushion_psf = P(1:4);
    P_cushion_Pa = P_cushion_psf' * psf_to_pa;

end

function F = residual_function(P, S, N_vec, dzdt_vec, Area)
    % P: [4x1] 气室压强 [P_FR, P_FL, P_RR, P_RL]
    % S: [4x1] 泄流面积
    % N_vec: [4x1] 风机转速
    
    % --- 1. 风机进气流量 (Q_FAN) ---
    % 修改逻辑：每个风机直接对着一个气室吹
    % 原公式中 term 依赖于 P_duct - 300。
    % 简化处理：我们将风机特性直接作用于气室压强。
    % 为了保持原有风机的量级，我们保留原公式的系数，但把压强换成当前气室压强。
    % 注意：如果压强过低，原公式 (-1280*sqrt) 可能会导致流量过大。
    % 这里假设风机出口直接就是气室。
    
    Q_FAN = zeros(4,1);
    for i = 1:4
        % 使用原有的风机特性曲线参数，但 P(5) 替换为 P(i)
        % 300 是原系统的一个参考背压，如果气室压强较低(100左右)，
        % P(i)-300 是负数，signed_sqrt 会处理它，产生正向流量。
        term1 = -1280 * signed_sqrt(P(i) - 300);
        term2 = -31.6 * (P(i) - 300);
        Q_FAN(i) = (term1 + term2) * N_vec(i) / 2000; 
    end

    % --- 2. 气室间横流 (Q_IC) ---
    % 拓扑结构：1(FR) -- 2(FL) -- 3(RR) -- 4(RL) -- 1(FR) (假设环形或依据原代码逻辑)
    % 原代码逻辑：
    % Q_IC1: 4 -> 1
    % Q_IC2: 1 -> 2
    % Q_IC3: 2 -> 3
    % Q_IC4: 3 -> 4
    
    Q_IC1 = 675 * signed_sqrt(P(4) - P(1));
    Q_IC2 = 338 * signed_sqrt(P(1) - P(2));
    Q_IC3 = 675 * signed_sqrt(P(2) - P(3));
    Q_IC4 = 338 * signed_sqrt(P(3) - P(4));

    % --- 3. 底部泄流 (Q_Leak) ---
    Q_Leak = -S .* 14.5 .* signed_sqrt(P);

    % --- 4. 活塞效应/泵吸 (Q_Pump) ---
    Q_PUMP = Area * dzdt_vec;

    % --- 5. 流量平衡方程 (Flow Balance) ---
    % 此时不需要 F(5) 和 F(6) 了，只需要 4 个方程
    F = zeros(4, 1);
    
    % 节点 1 (前右): 进=Fan1+IC1, 出=IC2+Leak+Pump
    F(1) = Q_FAN(1) + Q_IC1 - Q_IC2 + Q_Leak(1) - Q_PUMP(1);
    
    % 节点 2 (前左): 进=Fan2+IC2, 出=IC3+Leak+Pump
    F(2) = Q_FAN(2) + Q_IC2 - Q_IC3 + Q_Leak(2) - Q_PUMP(2);
    
    % 节点 3 (后左): 进=Fan3+IC3, 出=IC4+Leak+Pump
    F(3) = Q_FAN(3) + Q_IC3 - Q_IC4 + Q_Leak(3) - Q_PUMP(3);
    
    % 节点 4 (后右): 进=Fan4+IC4, 出=IC1+Leak+Pump (注意原代码 Q_IC1 是 4->1，所以这里出是 IC1 的反向? 
    % 不，原代码 F(4) = Q_INC4 + Q_IC4 - Q_IC1...
    % Q_IC4 是 3->4 (进入4), Q_IC1 是 4->1 (离开4)
    F(4) = Q_FAN(4) + Q_IC4 - Q_IC1 + Q_Leak(4) - Q_PUMP(4);

end

function J = compute_jacobian(P, S, N_vec, dzdt_vec, Area)
    eps_pert = 1e-4; % 扰动步长
    n = 4; % 现在是4维
    J = zeros(n, n);
    F0 = residual_function(P, S, N_vec, dzdt_vec, Area);
    
    for k = 1:n
        P_tmp = P;
        P_tmp(k) = P_tmp(k) + eps_pert;
        F_pert = residual_function(P_tmp, S, N_vec, dzdt_vec, Area);
        J(:, k) = (F_pert - F0) / eps_pert;
    end
end

%% ================= 辅助函数 =================
function val = signed_sqrt(x)
    % 带符号的开方函数
    val = sqrt(abs(x)) .* sign(x);
end
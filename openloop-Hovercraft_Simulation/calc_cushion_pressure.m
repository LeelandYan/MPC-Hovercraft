function [P_cushion_Pa, P_internal_state] = calc_cushion_pressure(N_FAN, S_leak_ft2, dzdt_ft_s, P_guess_psf)
% CALC_CUSHION_PRESSURE 计算气垫船气室压强
%
% 输入参数:
%   N_FAN       : [1x2] 向量, 左右风机转速 (RPM), 例如 [1800, 1800]
%   S_leak_ft2  : [1x4] 向量, 4个气室的底部泄流面积 (平方英尺)
%   dzdt_ft_s   : 标量, 船体垂直升沉速度 (英尺/秒) 
%   P_guess_psf : [6x1] 向量, 压强初值猜测 (psf)
%
% 输出参数:
%   P_cushion_Pa   : [1x4] 向量, 底部4个气室的最终压强 (帕斯卡 Pa)

%   P_internal_state : [6x1] 向量, 包含风机管道和气室的完整压强状态 (psf)

    %% 常数定义 
    Area_cushion_ft2 = 800; % 单个气室面积 (对应原文 74.32 m^2)
    psf_to_pa = 47.8803;    % 单位转换系数
    
    % 求解器配置
    max_iter = 30;
    tol = 20; % 容差
    lambda = 0.6; % 牛顿法步长因子

    %% 检查输入
    if nargin < 4 || isempty(P_guess_psf)
        P = ones(6,1) * 10; 
    else
        P = P_guess_psf;
    end
    
    N1 = N_FAN(1);
    N2 = N_FAN(2);
    
    % 泄流面积向量 S
    S = S_leak_ft2; 

    dzdt_vec = dzdt_ft_s(:);

    %% 牛顿迭代求解
    converged = false;
    
    for k = 1:max_iter
        % 计算残差 F
        F = residual_function(P, S, N1, N2, dzdt_vec, Area_cushion_ft2);
        
        % 检查收敛
        if norm(F) < tol
            converged = true;
            break; 
        end
        
        % 计算雅可比矩阵
        J = compute_jacobian(P, S, N1, N2, dzdt_vec, Area_cushion_ft2);
        
        % 更新步长
        delta = -J \ F;
        P = P + lambda * delta;
        
        P = max(P, 0.1); 
    end
    
    if ~converged
        warning('Pressure solver did not converge within max iterations.');
    end

    %% 输出处理
    P_internal_state = P; % 保持 psf 单位用于下一次迭代
    
    % 提取底部4个气室的压强并转换为 Pa
    P_cushion_psf = P(1:4);
    P_cushion_Pa = P_cushion_psf' * psf_to_pa;

end

function F = residual_function(P, S, N1, N2, dzdt_vec, Area)
    % P(1-4): 气室压强, P(5-6): 风机管道压强
    
    % 对应原文 Q_FAN 公式
    term1_1 = -1280 * signed_sqrt(P(5) - 300);
    term1_2 = -31.6 * (P(5) - 300);
    Q_FAN1 = (term1_1 + term1_2) * N1 / 2000;
    
    term2_1 = -1280 * signed_sqrt(P(6) - 300);
    term2_2 = -31.6 * (P(6) - 300);
    Q_FAN2 = (term2_1 + term2_2) * N2 / 2000;

    % 风机管道 -> 气室
    Q_INC1 = 589 * signed_sqrt(P(5) - P(1));
    Q_INC2 = 589 * signed_sqrt(P(5) - P(2));
    Q_INC3 = 589 * signed_sqrt(P(6) - P(3));
    Q_INC4 = 589 * signed_sqrt(P(6) - P(4));

    % 姿态控制喷管
    Q_NOZ1 = -346 * signed_sqrt(P(5));
    Q_NOZ2 = -346 * signed_sqrt(P(6));

    % 气室间横流
    Q_IC1 = 675 * signed_sqrt(P(4) - P(1));
    Q_IC2 = 338 * signed_sqrt(P(1) - P(2));
    Q_IC3 = 675 * signed_sqrt(P(2) - P(3));
    Q_IC4 = 338 * signed_sqrt(P(3) - P(4));

    % 底部泄流
    Q1 = -S(1) * 14.5 * signed_sqrt(P(1));
    Q2 = -S(2) * 14.5 * signed_sqrt(P(2));
    Q3 = -S(3) * 14.5 * signed_sqrt(P(3));
    Q4 = -S(4) * 14.5 * signed_sqrt(P(4));

    % 泵吸
    Q_PUMP = Area * dzdt_vec;

    % 流量平衡方程
    F = zeros(6, 1);
    
    % 风机管道节点平衡
    F(5) = -Q_INC1 - Q_INC2 + Q_NOZ1 + Q_FAN1;
    F(6) = -Q_INC3 - Q_INC4 + Q_NOZ2 + Q_FAN2;
    
    % 气室节点平衡
    F(1) = Q_INC1 + Q_IC1 - Q_IC2 + Q1 - Q_PUMP(1);
    F(2) = Q_INC2 + Q_IC2 - Q_IC3 + Q2 - Q_PUMP(2);
    F(3) = Q_INC3 + Q_IC3 - Q_IC4 + Q3 - Q_PUMP(3);
    F(4) = Q_INC4 + Q_IC4 - Q_IC1 + Q4 - Q_PUMP(4);
end

function J = compute_jacobian(P, S, N1, N2, dzdt_vec, Area)
    eps_pert = 1e-4; % 扰动步长
    n = 6;
    J = zeros(n, n);
    F0 = residual_function(P, S, N1, N2, dzdt_vec, Area);
    
    for k = 1:n
        P_tmp = P;
        P_tmp(k) = P_tmp(k) + eps_pert;
        F_pert = residual_function(P_tmp, S, N1, N2, dzdt_vec, Area);
        J(:, k) = (F_pert - F0) / eps_pert;
    end
end

%% ================= 辅助函数 =================
function val = signed_sqrt(x)
    % 带符号的开方函数，处理负压差的情况
    val = sqrt(abs(x)) * sign(x);
end
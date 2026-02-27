function [opt_rudder_deg, U_opt] = run_hovercraft_nmpc(X_curr, ref_data, env_data, prev_U, Np, Ts)
% 输入参数:
%   X_curr   : 当前状态向量 [x, y, psi, u, v, r]' 
%   ref_data : 导引参考数据 [y_e_curr, psi_ref, pi_h]' (横向误差, 期望航向, 路径切向角)
%   env_data : 外部环境与干扰 [phi, theta, wind_speed, wind_dir, prop_rpm]'
%   prev_U   : 上一控制周期的最优控制序列 
%   Np       : 预测步数 
%   Ts       : MPC 离散时间步长 
%
% 输出参数:
%   opt_rudder_deg : 当前时刻应下发的最佳舵角 (度)
%   U_opt          : 优化后的整个控制序列，留作下一次的热启动

    %% 1. 权重与物理约束定义
    % 代价函数权重矩阵
    Q_ye = 0;     % 横向误差惩罚
    Q_psie = 100000;  % 航向误差惩罚
    R_du = 1;      % 舵角变化率惩罚

    % 物理约束 
    rudder_max = 30;              % 最大绝对舵角 (度)
    rudder_rate_max_sec = 15;     % 最大打舵速率 (度/秒)
    
    % 转换到离散步长内的增量约束
    max_delta_u = rudder_rate_max_sec * Ts; 

    %% 2. 构建优化问题的边界与线性约束
    % 上下界约束 
    lb = -rudder_max * ones(Np, 1);
    ub =  rudder_max * ones(Np, 1);

    % 变化率约束转换矩阵 (A * U <= b)
    % U(i) - U(i-1) <= max_delta_u  =>  U(i) - U(i-1) <= dU
    % U(i-1) - U(i) <= max_delta_u  => -U(i) + U(i-1) <= dU
    A = zeros(2*Np, Np);
    b = max_delta_u * ones(2*Np, 1);
    

    u_last = prev_U(1); 
    A(1, 1) = 1;  b(1) = max_delta_u + u_last;
    A(2, 1) = -1; b(2) = max_delta_u - u_last;
    

    row = 3;
    for i = 2:Np
        A(row, i) = 1;   A(row, i-1) = -1;  b(row) = max_delta_u;
        row = row + 1;
        A(row, i) = -1;  A(row, i-1) = 1;   b(row) = max_delta_u;
        row = row + 1;
    end

    %% 3. 求解器配置 (使用 SQP 算法)
    options = optimoptions('fmincon', ...
        'Algorithm', 'sqp', ...          % 匹配文献中的算法
        'Display', 'none', ...           % 关闭内部打印以提升速度
        'MaxIterations', 50, ...         % 限制最大迭代次数
        'StepTolerance', 1e-3);          

    %% 4. 执行非线性优化
    % 目标函数句柄
    cost_func = @(U) compute_cost(U, X_curr, ref_data, env_data, Np, Ts, Q_ye, Q_psie, R_du, u_last);
    
    % 使用 fmincon 求解
    [U_opt, ~, exitflag] = fmincon(cost_func, prev_U, A, b, [], [], lb, ub, [], options);

    % 如果求解失败，保持上一次的舵角，保证系统安全
    if exitflag <= 0
        U_opt = prev_U; 
    end

    % 提取当前时刻最优控制量
    opt_rudder_deg = U_opt(1);
end

%% ================== 内部辅助函数：代价函数计算 ==================
function J = compute_cost(U, X, ref, env, Np, Ts, Q_ye, Q_psie, R_du, u_last)
    J = 0;
    x_k = X;
    
    % 解析参考数据与环境数据
    y_e_curr = ref(1);
    psi_ref  = ref(2);
    pi_h     = ref(3); % 路径切向角
    
    prev_u = u_last;

    % 在预测区间内滚动推演
    for i = 1:Np
        % 1. 提取当前预测步的试探控制量
        u_k = U(i);
        
        % 2. 计算控制增量惩罚 (减少机械磨损)
        delta_u = u_k - prev_u;
        J = J + R_du * (delta_u)^2;
        prev_u = u_k;
        

        dt_inner = 0.05; 
        num_inner = round(Ts / dt_inner);
        for inner = 1:num_inner
            % 3. 调用 3-DOF 预测模型计算状态导数
            dX = predict_3dof_model(x_k, u_k, env);
            
            % 4. 前向欧拉积分更新状态 (用 0.05s 的小步长推进)
            x_k = x_k + dX * dt_inner;
            
            % 提取临时预测状态
            psi_pred = x_k(3);
            u_pred   = x_k(4);
            v_pred   = x_k(5);
            
            % 6. 更新未来的横向误差 
            y_e_dot = u_pred * sin(psi_pred - pi_h) + v_pred * cos(psi_pred - pi_h);
            y_e_curr = y_e_curr + y_e_dot * dt_inner;
        end
        % ================================================================
        
        % 7. 计算航向误差 
        psi_err = angdiff(psi_ref, x_k(3));
        
        % 8. 累加预测步的代价值 
        J = J + Q_ye * (y_e_curr)^2 + Q_psie * (psi_err)^2;
    end

end

%% 
% 预测模型
function dX = predict_3dof_model(X, rudder_deg, env)
    % 状态解析
    % x_pos = X(1); y_pos = X(2); % 绝对位置在动力学内部不直接受力
    psi = X(3); u = X(4); v = X(5); r = X(6);
    
    % 环境与测量解析
    phi = env(1); theta = env(2); 
    wind_speed = env(3); wind_dir = env(4); prop_rpm = env(5);
    
    % 常数定义 (与你的 6-DOF 模型一致)
    FT2M = 0.3048; LBF2N = 4.44822; SLUG2KG = 14.5939; g_SI = 9.80665; rho_air_SI = 1.225;
    m_kg = 10879.5 * SLUG2KG;
    Izz = 2.057e7 * SLUG2KG * (FT2M^2);
    
    % 几何参数
    X_air_pos_si = 3; Z_air_pos_si = -3.0;
    X_rudder_pos_si = -67.1 * FT2M; Z_rudder_pos_si = -8.0 * FT2M;
    RSAREA_si = 47.2 * (FT2M^2); SCAREA_si = 3200 * (FT2M^2); FCAREA_si = 836 * (FT2M^2);
    DUCT_AREA_si = 123 * (FT2M^2);
    
    %% --- 空气动力计算 ---
    V_wind_north = wind_speed * cos(wind_dir + pi);
    V_wind_east  = wind_speed * sin(wind_dir + pi);
    u_wind_ship = V_wind_north * cos(psi) + V_wind_east * sin(psi);
    v_wind_ship = -V_wind_north * sin(psi) + V_wind_east * cos(psi);
    u_rel = u_wind_ship - u; v_rel = v_wind_ship - v;
    apparent_wind_velocity_si = sqrt(u_rel^2 + v_rel^2) + 0.001;
    abs_beta = abs(rad2deg(atan2(-v_rel, -u_rel)));
    
    % 阻力系数
    if abs_beta <= 180, SDRAG = 0.01385 * abs_beta - 7.69e-5 * abs_beta^2; else, SDRAG = 0; end
    if v_rel > 0, SDRAG = -SDRAG; end % 修正相对风向的正负逻辑
    
    if abs_beta < 40
        FDRAG = -2.22e-4 * abs_beta^2 + 3.33e-3 * abs_beta + 0.5;
    elseif abs_beta < 88
        FDRAG = -2.48e-5 * abs_beta^3 + 5e-3 * abs_beta^2 - 0.324 * abs_beta + 5.835;
    else
        FDRAG = 1.04e-4 * abs_beta^2 - 3.518e-2 * abs_beta + 2.39;
    end
    
    q_bar = 0.5 * rho_air_SI * apparent_wind_velocity_si^2;
    X_air_si = -FDRAG * FCAREA_si * q_bar;
    Y_air_si = -SDRAG * SCAREA_si * q_bar;
    Mz_air_si = Y_air_si * X_air_pos_si;

    %% --- 推力与舵力计算 ---
    XWIND_val = -u_rel / FT2M;
    prop_angle = 15; % 根据原模型预设
    Thrust_base = (338 * prop_angle + 4.36 * prop_angle^2);
    Thrust_loss = (1.43 * prop_angle * XWIND_val + 0.1715 * XWIND_val^2);
    Thrust_one_si = (Thrust_base - Thrust_loss) * (prop_rpm/1250)^2 * LBF2N;
    X_prop_si = 2 * Thrust_one_si;
    
    VDUCT2_si = (apparent_wind_velocity_si * cos(atan2(-v_rel, -u_rel)))^2 + Thrust_one_si / (rho_air_SI * DUCT_AREA_si);
    if VDUCT2_si < 0, VDUCT2_si = 0; end
    PDUCT_si = 0.5 * rho_air_SI * VDUCT2_si;
    
    CLIFT_R = 0.053 * rudder_deg;
    CDRAG_R = 0.422e-3 * rudder_deg^2;
    X_rudder_si = -CDRAG_R * RSAREA_si * 2 * PDUCT_si;
    Y_rudder_si = CLIFT_R * RSAREA_si * 2 * PDUCT_si;
    Mz_rudder_si = X_rudder_pos_si * Y_rudder_si;
    
    %% --- 水动力阻力 ---
    Rho_water = 1025;
    X_skirt_si = -0.5 * Rho_water * 0.25 * (50 * FT2M^2) * u * abs(u);
    Y_skirt_si = -0.5 * Rho_water * 3 * (80 * FT2M^2) * v * abs(v);
    Mz_skirt_si = - (2.77e6 * LBF2N * FT2M) * r;

    %% --- 引入实时测量的重力分量 (降维补偿核心) ---
    X_G_si = -m_kg * g_SI * sin(theta);
    Y_G_si =  m_kg * g_SI * sin(phi);

    %% --- 运动学与动力学积分 ---
    dx = u * cos(psi) - v * sin(psi);
    dy = u * sin(psi) + v * cos(psi);
    dpsi = r;
    
    Fx = X_rudder_si + X_air_si + X_G_si + X_prop_si + X_skirt_si;
    Fy = Y_rudder_si + Y_air_si + Y_G_si + Y_skirt_si;
    Mz = Mz_air_si + Mz_rudder_si + Mz_skirt_si;
    
    udot = Fx/m_kg + r*v;
    vdot = Fy/m_kg - r*u;
    rdot = Mz/Izz;

    dX = [dx; dy; dpsi; udot; vdot; rdot];
end

%% ================== 内部辅助函数：计算最小角度差 ==================
function d = angdiff(th1, th2)
    % 将角度差转换到 [-pi, pi] 之间，防止 359度 和 1度 算成 358度 的巨大误差
    d = th1 - th2;
    d = mod(d + pi, 2*pi) - pi;
end
function [dXdt, P_cushion_out] = model_jeff_b(t, X, cmd_rudder_angle, wind_param) 
    %% --- 单位转换系数 ---
    FT2M = 0.3048;          % feet to meters
    SLUG2KG = 14.5939;      % slugs to kg
    LBF2N = 4.44822;        % pound-force to Newtons
    rho_air_SI = 1.225;     % kg/m^3
    g_SI = 9.80665;         % m/s^2
  
    M2FT = 1/FT2M;

    %% --- 定义状态量 ---
    phi   = X(4); % Roll (横摇)
    theta = X(5); % Pitch (纵摇)
    psi   = X(6); % Yaw (艏向)
    u = X(7);  % 纵向速度
    v = X(8);  % 横向速度
    w = X(9);  % 垂向速度
    p = X(10); % 横摇角速度
    q = X(11); % 纵摇角速度
    r = X(12); % 艏摇角速度
    
    %% --- 气垫船参数 ---
    m_slugs = 10879.5;
    Ixx_imp = 5.672e6;  
    Iyy_imp = 1.629e7;
    Izz_imp = 2.057e7;
    
    % 转换为国际单位制
    m_kg = m_slugs * SLUG2KG;
    Ixx = Ixx_imp * SLUG2KG * (FT2M^2);
    Iyy = Iyy_imp * SLUG2KG * (FT2M^2);
    Izz = Izz_imp * SLUG2KG * (FT2M^2);
    
    % 空气阻力作用点
    X_air_pos_si = 3;    % 纵向偏移 (x)
    Z_air_pos_si = -3.0 ;    % 垂向偏移 (z)

    
    % 几何参数
    X_rudder_pos_si = -67.1 * FT2M;       % 舵 X 位置
    Z_rudder_pos_si = -8.0 * FT2M;        % 舵 Z 位置
    RSAREA_si = 47.2 * (FT2M^2);  
    SCAREA_si = 3200 * (FT2M^2);  
    FCAREA_si = 836 * (FT2M^2);   
    DUCT_AREA_si = 123 * (FT2M^2);
    Z_propeller_pos_si = -3; % 螺旋桨安装高度


    %% 控制输入 
    propeller_angle = 15;

%     propeller_speed_rpm = 1200;
    if t > 10
%         propeller_speed_rpm = 1200;
        propeller_speed_rpm = 0;
    else
        propeller_speed_rpm = 0;
    end

    
    if t > 100
        rudder_angle_deg = cmd_rudder_angle;
    else
        rudder_angle_deg = 0;
    end

    rudder_angle_rad = deg2rad(rudder_angle_deg);
    
%     true_wind_speed_si = 0;       % 真风速
%     true_wind_direction_si = deg2rad(90);   % 真风来向
    [true_wind_speed_si, true_wind_direction_si] = Get_Wind_Step(t, wind_param);



    %% 气垫几何参数
    % 气垫分布几何 
    L_cush = 38.5 * FT2M; % 单气室长11.73m
    W_cush = 17.5 * FT2M; % 单气室宽5.3m

    % 气室中心坐标
    Pos_cush = [
         L_cush/2,  W_cush/2; 
         L_cush/2, -W_cush/2; 
        -L_cush/2,  W_cush/2; 
        -L_cush/2, -W_cush/2
    ];
    
    % 围裙与风机参数
    L_per_cushion_ft = 140;         % 单个气室围裙下缘长度
    Area_cushion_ft2 = 800;         % 单个气室面积
    SD_equilibrium_ft = 4.5;        % 平衡围裙深度
    P_equilibrium_psf = 109;        % 平衡压强
    H_base_m = 5.0 * FT2M;          % 硬结构离水面基准高度
    
    C_SKRT = 0.0112;                % 围裙刚度系数
    TC = 8.0;                       % 围裙响应时间常数

    Target_RPM_1 = 1500;
    Run_Up_Time = 60.0;
    if t < Run_Up_Time
        Current_RPM_1 = Target_RPM_1 * (1 - exp(-0.5 * t/1));  % 指数上升
    else
        Current_RPM_1 = Target_RPM_1;
    end

    Target_RPM_2 = 2000;
    Run_Up_Time = 60.0;
    if t < Run_Up_Time
        Current_RPM_2 = Target_RPM_2 * (1 - exp(-0.5 * t/1));  % 指数上升
    else
        Current_RPM_2 = Target_RPM_2;
    end
    
   
    N_FAN_FL = Current_RPM_2; % 前左
    N_FAN_FR = 800;% Current_RPM_1; % 前右
    N_FAN_RL = 800;% Current_RPM_1; % 后左
    N_FAN_RR = 800;% Current_RPM_1; % 后右
    Fan_RPMs = [N_FAN_FR, N_FAN_FL, N_FAN_RR, N_FAN_RL];
    

%% 空气阻力计算    
    % --- 相对风计算(Apparent Wind) Eq.60-63 ---
    % 将真风速度分解到固定坐标系 (NED)
    V_wind_north = true_wind_speed_si * cos(true_wind_direction_si + pi);
    V_wind_east  = true_wind_speed_si * sin(true_wind_direction_si + pi);
    
    % 将真风速度分量转换到船体坐标系
    u_wind_ship = V_wind_north * cos(psi) + V_wind_east * sin(psi);
    v_wind_ship = -V_wind_north * sin(psi) + V_wind_east * cos(psi);
    
    % 计算相对风速
    % 相对风 = 真风分量 - 船速
    u_rel = u_wind_ship - u;
    v_rel = v_wind_ship - v;
    
    % 相对风速
    apparent_wind_velocity_si = sqrt(u_rel^2 + v_rel^2) + 0.001; % 防止除零
    
    % 计算相对风角 (Apparent Wind Angle, Beta)
    apparent_wind_angle_rad = atan2(-v_rel, -u_rel);  % 相对风来向角（弧度）：从船首起算，顺时针为正
    apparent_wind_angle_deg = rad2deg(apparent_wind_angle_rad);
    abs_beta  = abs(apparent_wind_angle_deg);
    
    % --- 气动系数计算 Eq.70 ---
    % 侧向力系数 Eq.70 
    if abs_beta <= 180
        SDRAG = 0.01385 * abs_beta - 7.69e-5 * abs_beta^2;
    else
        SDRAG = 0;
    end
   
    if apparent_wind_angle_deg < 0, SDRAG = -SDRAG; end
    
    % 正面阻力系数 Eq.71 
    if abs_beta < 40
        FDRAG = -2.22e-4 * abs_beta^2 + 3.33e-3 * abs_beta + 0.5;
    elseif abs_beta < 88
        FDRAG = -2.48e-5 * abs_beta^3 + 5e-3 * abs_beta^2 ...
                      - 0.324 * abs_beta + 5.835;
    else
        FDRAG = 1.04e-4 * abs_beta^2 - 3.518e-2 * abs_beta + 2.39;
    end
    
    % 艏摇力矩系数 Eq.72 
    if abs_beta < 60
        YDRAG = 1.67e-3 * abs_beta;
    elseif abs_beta <= 120
        YDRAG = 5.17e-7 * abs_beta^3 - 1.23e-4 * abs_beta^2 ...
                      + 5.82e-3 * abs_beta + 8.27e-2;
    else
        YDRAG = 1.67e-3 ; % 保持常数
    end

    % --- 空气阻力与力矩计算 Eq.73-77 ---
    q_bar = 0.5 * rho_air_SI * apparent_wind_velocity_si^2; % 动压

    % 空气阻力计算
    % Eq.73: XBDRAG (纵向空气阻力) 
    X_air_si = -FDRAG * FCAREA_si * q_bar;
    
    % Eq.74: YBDRAG (横向空气阻力) 
    Y_air_si = -SDRAG * SCAREA_si * q_bar; 

    % 力矩 (Moments)
    % Eq.75: Pitch Moment
    My_air_si = X_air_si * Z_air_pos_si;  
    
    % Eq.76: Roll Moment 
    Mx_air_si  = -Y_air_si * Z_air_pos_si;
    
    % Eq.77: Yaw Moment 
%     Mz_air_si = (YDRAG * SCAREA_si * LCUSH_si * q_bar) ...
%              + (Y_air_si * X_air_pos_si);
     Mz_air_si = Y_air_si * X_air_pos_si;
    

    %% 推进与舵力计算
    XWIND_ft_s = -u_rel * M2FT; 
    XWIND_val = XWIND_ft_s; 
    Thrust_base = (338 * propeller_angle + 4.36 * propeller_angle^2);
    Thrust_loss = (1.43 * propeller_angle * XWIND_val + 0.1715 * XWIND_val^2);
    Thrust_one_lbf = (Thrust_base - Thrust_loss) * (propeller_speed_rpm/1250)^2;
    
    % 转换为牛顿
    Thrust_one_si = Thrust_one_lbf * LBF2N;
    
    % 总推力
    X_propeller_thrust_si = 2 * Thrust_one_si;

    My_propeller_si = X_propeller_thrust_si * Z_propeller_pos_si;
    
    % 舵力
    VDUCT2_si = (apparent_wind_velocity_si * cos(apparent_wind_angle_rad))^2 + Thrust_one_si / (rho_air_SI * DUCT_AREA_si);
    if VDUCT2_si < 0, VDUCT2_si = 0; end
    PDUCT_si = 0.5 * rho_air_SI * VDUCT2_si;
    
    CLIFT_R = 0.053 * rudder_angle_deg;
    CDRAG_R = 0.422e-3 * rudder_angle_deg^2;
    
    X_rudder_si = -CDRAG_R * RSAREA_si * 2 * PDUCT_si;
    Y_rudder_si = CLIFT_R * RSAREA_si * 2 * PDUCT_si;
    
    % 舵力矩
    % Mx
    Mx_rudder_si  = -Z_rudder_pos_si * Y_rudder_si;
    % My
    My_rudder_si  = Z_rudder_pos_si * X_rudder_si;
    % Mz
    Mz_rudder_si  = X_rudder_pos_si * Y_rudder_si;
    
%% 气垫力 
    % 记忆上一时刻的围裙状态和压强
    persistent P_last_psf SD_curr_ft t_last
    
    % 初始化 
    if isempty(P_last_psf) || t == 0
        P_last_psf = zeros(4,1); 
        SD_curr_ft = ones(4,1) * SD_equilibrium_ft; 
        t_last = t;
    end
    
    % 计算动态时间步长 dt
    dt_step = t - t_last;
    if dt_step < 0, dt_step = 0; end
    if dt_step > 0.1, dt_step = 0.01; end % 限制最大步长
    t_last = t;

    % 计算各气室物理高度
    z_ned = X(3); % 垂向位移 (向下为正)
    h_hull_ft = zeros(4,1);
    
    % 【修改点1】初始化一个向量，用来存4个气室的局部垂向速度
    dzdt_vec_ft = zeros(4,1); 
    
    for i = 1:4
        % 获取气室中心相对于重心的坐标 (力臂)
        x_arm = Pos_cush(i,1);
        y_arm = Pos_cush(i,2);

        % 计算高度
        % h = 基准 - 沉降 + 纵摇影响 - 横摇影响
        h_local_m = (H_base_m - z_ned) + x_arm*sin(theta) - y_arm*sin(phi);
        h_hull_ft(i) = h_local_m * M2FT; % 转为 ft


   
        dzdt_local_si = -w + q * x_arm - p * y_arm; 
        dzdt_vec_ft(i) = dzdt_local_si * M2FT; 

        % 围裙
        % 计算压强偏差 -> 目标深度
        Y_n = P_equilibrium_psf - P_last_psf(i);
        Y_n = max(min(Y_n, 66.9), -66.9); % 限幅
        SDPRO = 4.5 + C_SKRT * Y_n - 8.325e-7 * (Y_n^3);
        
        % 状态积分: 更新围裙长度
        SDDT = (SDPRO - SD_curr_ft(i)) / TC;
        SD_curr_ft(i) = SD_curr_ft(i) + SDDT * dt_step;
        
        % 计算物理气隙 (Gap)
        SAW = h_hull_ft(i) - SD_curr_ft(i);
        CLR = max(SAW, 1e-4); % 最小气隙保护
        
        % 泄流面积
        S_vec_ft2(i) = L_per_cushion_ft * CLR;
    end


    [P_cushion_Pa, P_internal_psf] = calc_cushion_pressure(...
        Fan_RPMs, S_vec_ft2', dzdt_vec_ft, P_last_psf);
    
    % 更新记忆变量
    P_last_psf = P_internal_psf;

    % 气垫力与力矩
    % 气垫升力
    F_lift_N = P_cushion_Pa * (Area_cushion_ft2 * FT2M^2); 
    
    % 垂向力(向上为负)
    Z_cursion_si = -sum(F_lift_N); 
    
    % 力矩 
    My_cursion_si = 0; % Pitch 
    Mx_cursion_si = 0; % Roll 
    
    for i = 1:4
        x_arm = Pos_cush(i,1);
        y_arm = Pos_cush(i,2);
        F_i = F_lift_N(i);
        
        % 纵摇力矩
        My_cursion_si = My_cursion_si + (F_i * x_arm);
        % 横摇力矩 
        Mx_cursion_si = Mx_cursion_si - (F_i * y_arm);
    end


    %% 航行面阻力
    CD_skirt = 0.25;
    Area_wet_surge = 50 * (FT2M^2);
    Area_wet_sway  = 80 * (FT2M^2);
    Rho_water = 1025;
    
    X_skirt_si = -0.5 * Rho_water * 0.25 * Area_wet_surge * u * abs(u);
    Y_skirt_si = -0.5 * Rho_water * 3 * Area_wet_sway  * v * abs(v);


    Z_arm_skirt = 1.5; % 阻力作用点到重心的垂直距离 (米)
    Mx_skirt_moment = -Z_arm_skirt * Y_skirt_si;
    Roll_Damping_Coeff = 5e5; 
    Mx_damping = -Roll_Damping_Coeff * p;

    Pitch_Damping_Coeff = 2e6; 
    My_damping = - Pitch_Damping_Coeff * q;
    Z_arm_effective = 1.0;
    My_skirt_moment = Z_arm_effective * X_skirt_si;
    
    YAWDC_coeff_si = 2.77e6 * LBF2N * FT2M;
    Mz_skirt_si = -YAWDC_coeff_si * r;    
    
    %% --- 7. 总力与力矩汇总 ---
    % 重力在各轴分量，近似线性
    X_G_si = -m_kg * g_SI * theta;
    Y_G_si =  m_kg * g_SI * phi;
    Z_G_si =  m_kg * g_SI;
    
    % 合力与合力矩
    Fx = X_rudder_si + X_air_si  + X_G_si + X_propeller_thrust_si + X_skirt_si;
    Fy = Y_rudder_si + Y_air_si  + Y_G_si + Y_skirt_si;
    Fz = Z_cursion_si + Z_G_si;
    
    Mx = Mx_cursion_si + Mx_air_si + Mx_rudder_si + Mx_skirt_moment + Mx_damping;
    My = My_cursion_si + My_air_si + My_rudder_si + My_propeller_si + My_skirt_moment + My_damping;
    Mz = Mz_air_si + Mz_rudder_si + Mz_skirt_si;
    
    %% --- 8. 动力学方程 ---    
    udot = Fx/m_kg - q*w + r*v;
    vdot = Fy/m_kg - r*u + p*w;
    wdot = Fz/m_kg - p*v + q*u;
    
    pdot = (Mx - (Izz - Iyy)*q*r) / Ixx;
    qdot = (My - (Ixx - Izz)*r*p) / Iyy;
    rdot = (Mz - (Iyy - Ixx)*p*q) / Izz;
    
    %% --- 9. 运动学方程 ---
    dx = u * cos(psi) - v * sin(psi);
    dy = u * sin(psi) + v * cos(psi);
    dz = w;
    
    dphi   = p;
    dtheta = q;
    dpsi   = r;
    
    %% --- 输出 ---
    dXdt = [dx; dy; dz; dphi; dtheta; dpsi; ...
            udot; vdot; wdot; pdot; qdot; rdot];

    P_cushion_out = P_cushion_Pa;

end
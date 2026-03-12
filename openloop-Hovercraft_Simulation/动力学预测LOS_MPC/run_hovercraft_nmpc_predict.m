function [opt_rudder_deg, U_opt] = run_hovercraft_nmpc_predict(X_curr, ref_data, env_data, prev_U, Np, Ts)
    import casadi.*
    
    persistent opti X U params solver_ready

    if isempty(solver_ready)
        opti = casadi.Opti();

        % 声明优化变量
        X = opti.variable(6, Np+1); 
        U = opti.variable(1, Np);

        % 声明实时变化的“参数”
        params.X0     = opti.parameter(6, 1); 
        params.ref    = opti.parameter(3, 1); 
        params.env    = opti.parameter(5, 1); 
        params.u_last = opti.parameter(1, 1); 

        y_e_curr_val = params.ref(1);
        pi_h_val     = params.ref(2);
        y_e_int_val  = params.ref(3);         % 接收积分项

        % 权重矩阵 
        Q_psi = 10000;      
        Q_r   = 10;      
        Q_u   = 10;         
        u_ref = 3.0;   
        
        R_du = 50;          
        R_u  = 1;           
        
        rudder_max = 30;              
        rudder_rate_max_sec =1;     
        max_delta_u = rudder_rate_max_sec * Ts; 

        % ILOS参数
        Delta_h = 150;      % 前视距离
        K_int   = 0.01;     % 积分增益系数，控制侧滑补偿的速度

        
        f_dynamics = get_casadi_dynamics();
        opti.subject_to(X(:, 1) == params.X0); 

        J = 0;
        dt_inner = 0.05; 
        num_inner = round(Ts / dt_inner);

        for i = 1:Np
            % 控制量约束 
            opti.subject_to(-rudder_max <= U(i) <= rudder_max);
            if i == 1
                delta_u = U(i) - params.u_last;
            else
                delta_u = U(i) - U(i-1);
            end
            opti.subject_to(-max_delta_u <= delta_u <= max_delta_u);
            
            % 累加控制指令代价
            J = J + R_du * (delta_u)^2 + R_u * (U(i))^2;

            % 状态预测
            x_temp = X(:, i);
            for inner = 1:num_inner
                dX = f_dynamics(x_temp, U(i), params.env);
                x_temp = x_temp + dX * dt_inner;
            end
            opti.subject_to(X(:, i+1) == x_temp);
            
            % 提取预测状态终端值
            x_pred_end   = X(1, i+1);
            y_pred_end   = X(2, i+1);
            psi_pred_end = X(3, i+1);
            u_pred_end   = X(4, i+1);
            r_pred_end   = X(6, i+1);
            
            % 计算预测的横向误差 y_e
            y_e_pred = y_e_curr_val - (x_pred_end - params.X0(1))*sin(pi_h_val) ...
                                    + (y_pred_end - params.X0(2))*cos(pi_h_val);
            

            psi_c_pred = pi_h_val - atan((y_e_pred + K_int * y_e_int_val) / Delta_h);
            
            % 计算追踪误差
            psi_err_cost = 2 * (1 - cos(psi_c_pred - psi_pred_end));
            u_err = u_ref - u_pred_end;
            
            J = J + Q_psi * psi_err_cost + Q_r * (r_pred_end)^2 + Q_u * (u_err)^2;
        end

        opti.minimize(J);

        p_opts = struct('expand', true, 'print_time', false);
        s_opts = struct('max_iter', 50, 'print_level', 0, 'acceptable_tol', 1e-3);
        opti.solver('ipopt', p_opts, s_opts);
        
        solver_ready = true;
    end

    opti.set_value(params.X0, X_curr);
    opti.set_value(params.ref, ref_data);
    opti.set_value(params.env, env_data);
    opti.set_value(params.u_last, prev_U(1));
    opti.set_initial(U, prev_U');
    
    X_init_guess = repmat(X_curr, 1, Np+1);
    opti.set_initial(X, X_init_guess);

    try
        sol = opti.solve();
        U_opt = sol.value(U)';
    catch
        U_opt = opti.debug.value(U)';
        if isnan(U_opt(1)), U_opt = prev_U; end
    end

    opt_rudder_deg = U_opt(1);
end

%% 辅助函数
function f = get_casadi_dynamics()
    import casadi.*
    
    x_sym = MX.sym('x_sym', 6); 
    u_sym = MX.sym('u_sym', 1); 
    env_sym = MX.sym('env_sym', 5); 
    
    psi = x_sym(3); u = x_sym(4); v = x_sym(5); r = x_sym(6);
    phi = env_sym(1); theta = env_sym(2); 
    wind_speed = env_sym(3); wind_dir = env_sym(4); prop_rpm = env_sym(5);

    FT2M = 0.3048; LBF2N = 4.44822; SLUG2KG = 14.5939; g_SI = 9.80665; rho_air_SI = 1.225;
    m_kg = 10879.5 * SLUG2KG;
    Izz = 2.057e7 * SLUG2KG * (FT2M^2);
    
    X_air_pos_si = 3; Z_air_pos_si = -3.0;
    X_rudder_pos_si = -67.1 * FT2M; Z_rudder_pos_si = -8.0 * FT2M;
    RSAREA_si = 47.2 * (FT2M^2); SCAREA_si = 3200 * (FT2M^2); FCAREA_si = 836 * (FT2M^2);
    DUCT_AREA_si = 123 * (FT2M^2);

    V_wind_north = wind_speed * cos(wind_dir + pi);
    V_wind_east  = wind_speed * sin(wind_dir + pi);
    u_wind_ship = V_wind_north * cos(psi) + V_wind_east * sin(psi);
    v_wind_ship = -V_wind_north * sin(psi) + V_wind_east * cos(psi);
    u_rel = u_wind_ship - u; 
    v_rel = v_wind_ship - v;
    apparent_wind_velocity_si = sqrt(u_rel^2 + v_rel^2 + 1e-6);
    
    epsilon = 1e-6; 
    abs_beta_rad = sqrt(atan2(-v_rel, -u_rel + 1e-6)^2 + epsilon);
    abs_beta = abs_beta_rad * (180/pi);
    
    SDRAG_val = 0.01385 * abs_beta - 7.69e-5 * abs_beta^2;
    smooth_sign_v = v_rel / sqrt(v_rel^2 + epsilon);
    SDRAG = if_else(abs_beta <= 180, SDRAG_val * (-smooth_sign_v), 0);
    
    FDRAG_1 = -2.22e-4 * abs_beta^2 + 3.33e-3 * abs_beta + 0.5;
    FDRAG_2 = -2.48e-5 * abs_beta^3 + 5e-3 * abs_beta^2 - 0.324 * abs_beta + 5.835;
    FDRAG_3 = 1.04e-4 * abs_beta^2 - 3.518e-2 * abs_beta + 2.39;
    FDRAG = if_else(abs_beta < 40, FDRAG_1, if_else(abs_beta < 88, FDRAG_2, FDRAG_3));
    
    q_bar = 0.5 * rho_air_SI * apparent_wind_velocity_si^2;
    X_air_si = -FDRAG * FCAREA_si * q_bar;
    Y_air_si = -SDRAG * SCAREA_si * q_bar;
    Mz_air_si = Y_air_si * X_air_pos_si;

    XWIND_val = -u_rel / FT2M;
    prop_angle = 15; 
    Thrust_base = (338 * prop_angle + 4.36 * prop_angle^2);
    Thrust_loss = (1.43 * prop_angle * XWIND_val + 0.1715 * XWIND_val^2);
    Thrust_one_si = (Thrust_base - Thrust_loss) * (prop_rpm/1250)^2 * LBF2N;
    X_prop_si = 2 * Thrust_one_si;
    
    VDUCT2_si_val = (-u_rel)^2 + Thrust_one_si / (rho_air_SI * DUCT_AREA_si);
    VDUCT2_si = 0.5 * (VDUCT2_si_val + sqrt(VDUCT2_si_val^2 + epsilon));
    PDUCT_si = 0.5 * rho_air_SI * VDUCT2_si;
    
    CLIFT_R = 0.053 * u_sym;
    CDRAG_R = 0.422e-3 * u_sym^2;
    X_rudder_si = -CDRAG_R * RSAREA_si * 2 * PDUCT_si;
    Y_rudder_si = CLIFT_R * RSAREA_si * 2 * PDUCT_si;
    Mz_rudder_si = X_rudder_pos_si * Y_rudder_si;

    Rho_water = 1025;
    X_skirt_si = -0.5 * Rho_water * 0.25 * (50 * FT2M^2) * u * sqrt(u^2 + epsilon);
    Y_skirt_si = -0.5 * Rho_water * 3 * (80 * FT2M^2) * v * sqrt(v^2 + epsilon);
    Mz_skirt_si = - (2.77e6 * LBF2N * FT2M) * r;

    X_G_si = -m_kg * g_SI * sin(theta); 
    Y_G_si =  m_kg * g_SI * sin(phi);

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
    
    f = Function('f_dynamics', {x_sym, u_sym, env_sym}, {dX});
end
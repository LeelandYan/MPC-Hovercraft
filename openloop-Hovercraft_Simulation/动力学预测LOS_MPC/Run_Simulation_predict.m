%% 气垫船6DOF仿真 
clear all; clc; close all;

%% 初始状态定义
% 状态向量 X = [x, y, z, phi, theta, psi, u, v, w, p, q, r]
x0 = -100;      % 北向位置 (m)
y0 = 0;      % 东向位置 (m)
z0 = 1.5;   % 垂向位移 (m)
phi0 = 0;    % 横摇角 (rad)
theta0 = 0;  % 纵摇角 (rad)
psi0 = deg2rad(1); % 初始艏向 (rad)，0度为正北

% 速度
u0 = 0.5;      % 纵向速度 (m/s)
v0 = 0;      % 横向速度 (m/s)
w0 = 0;      % 垂向速度 (m/s)
p0 = 0;      % 横摇角速度 (rad/s)
q0 = 0;      % 纵摇角速度 (rad/s)
r0 = 0;      % 艏摇角速度 (rad/s)

% 组装初始状态向量
X0 = [x0, y0, z0, phi0, theta0, psi0, u0, v0, w0, p0, q0, r0];

% 定义时间步长和总时长
dt = 0.05;              
T_end = 400;            % 仿真结束时间
t = 0:dt:T_end;         % 生成时间向量
num_steps = length(t);  % 总步数

%%
% 设定平均风速，风向
V_wind_mean = 5; 
Psi_wind_mean = deg2rad(180); 
wind_data = Init_Wind(V_wind_mean, Psi_wind_mean);

% 控制指令
rudder_angle = 0; % 舵角指令

%% LOS与MPC初始化 
% 航线 (单位米)
wpt.pos.x = [0, 1000, 1000, 0]';    % 北向 (North)
wpt.pos.y = [0, 1000, 2000, 3000]'; % 东向 (East)
% wpt.pos.x = [0, 500, 500]';    % 北向 (North)
% wpt.pos.y = [0, 500, 2000]'; % 东向 (East)
% wpt.pos.x = [0, 1500]';    % 北向 (North)
% wpt.pos.y = [0, 1500]'; % 东向 (East)

% LOS 参数设定
Delta_h = 180;      % 前视距离 (m)
R_switch = 150;    % 航点切换半径 (m)
clear LOSpsi;    

% MPC 参数设定 (修正非最小相位特性)
Ts_mpc = 1.0;      % 【修改】: 增大 MPC 控制周期至 1.0s
Np = 15;           % 【修改】: 增大预测步数至 15 步 (总视野 15s)
control_steps = round(Ts_mpc / dt); % 计算间隔步数 (1.0 / 0.05 = 20)
prev_U = zeros(Np, 1); % 初始化 MPC 控制序列
y_e_int = 0;

% 记录控制数据的数组
MPC_rudder_history = zeros(num_steps, 1);
MPC_rudder_history(1) = rudder_angle; % 记录第1步的舵角


%% 使用4阶龙格库塔法
% 状态矩阵 sol: [行=时间步, 列=12个状态量]
sol = zeros(num_steps, 12); 
sol(1, :) = X0;         % 填入初始状态

% 压强矩阵 P_hist_Pa: [行=时间步, 列=4个气室]
P_hist_Pa = zeros(num_steps, 4);

% 计算初始时刻的压强
[~, P_init] = model_jeff_b(t(1), X0', rudder_angle, wind_data); 
P_hist_Pa(1, :) = P_init(:)';

fprintf('进行气垫船6自由度仿真...\n');

% RK4循环
for k = 1 : num_steps - 1
    % 当前时刻和状态
    t_curr = t(k);
    X_curr = sol(k, :)'; % 取出为列向量 (12x1)

    if mod(k, round(10 / dt)) == 0
        fprintf('当前仿真进度: %5.1f / %.1f 秒 (已完成 %5.1f%%)\n', t_curr, T_end, (t_curr/T_end)*100);
    end
    
    % MPC与LOS
    if mod(k-1, control_steps) == 0
        % 提取当前观测状态
        x = X_curr(1); y = X_curr(2); 
        phi = X_curr(4); theta = X_curr(5); psi = X_curr(6);
        u = X_curr(7); v = X_curr(8); r = X_curr(12);
        
        % 1. 获取当前的 y_e 和 pi_h。
        [~, y_e, pi_h, ~] = LOSchi(x, y, Delta_h, R_switch, wpt);

        % 对横向误差进行数值积分
        y_e_int = y_e_int + y_e * Ts_mpc;

        % 2. 打包状态量和参考数据给MPC 
        X_mpc = [x; y; psi; u; v; r];
        ref_data = [y_e; pi_h; y_e_int];
        
        % 3. 获取当前环境数据
        [wind_speed, wind_dir] = Get_Wind_Step(t_curr, wind_data);
        prop_rpm = 1000; 
        env_data = [phi; theta; wind_speed; wind_dir; prop_rpm];
        
        % 4. 运行MPC求解器
        [opt_rudder_deg, prev_U] = run_hovercraft_nmpc_predict(X_mpc, ref_data, env_data, prev_U, Np, Ts_mpc);
        
        % 5. 下发舵角指令
        rudder_angle = opt_rudder_deg;
    end
   
    % RK4迭代
    [dX1, ~] = model_jeff_b(t_curr,          X_curr,              rudder_angle, wind_data);
    [dX2, ~] = model_jeff_b(t_curr + 0.5*dt, X_curr + 0.5*dt*dX1, rudder_angle, wind_data);
    [dX3, ~] = model_jeff_b(t_curr + 0.5*dt, X_curr + 0.5*dt*dX2, rudder_angle, wind_data);
    [dX4, ~] = model_jeff_b(t_curr + dt,     X_curr + dt*dX3,     rudder_angle, wind_data);
    X_next = X_curr + (dt / 6) * (dX1 + 2*dX2 + 2*dX3 + dX4);
    
    sol(k+1, :) = X_next';
    
    [~, P_out] = model_jeff_b(t(k+1), X_next, rudder_angle, wind_data);
    P_hist_Pa(k+1, :) = P_out(:)';

    % 记录舵角
    MPC_rudder_history(k+1) = rudder_angle;
end

fprintf('仿真完成。\n');

%% 数据解包 
% 位置 (m)
x_pos = sol(:,1);
y_pos = sol(:,2);
z_pos = sol(:,3); 

% 姿态 (弧度 -> 角度)
phi_deg   = rad2deg(sol(:,4));
theta_deg = rad2deg(sol(:,5));
psi_deg   = rad2deg(sol(:,6)); 

% 速度 (m/s)
u = sol(:,7);
v = sol(:,8);
w = sol(:,9);

% 角速度 (rad/s -> deg/s)
r_rate_deg = rad2deg(sol(:,12));


% 压强数据
P_FL = P_hist_Pa(:,1); % 前左  
P_FR = P_hist_Pa(:,2); % 前右  
P_RR = P_hist_Pa(:,3); %      
P_RL = P_hist_Pa(:,4); %      

P_static_si = 5218; % 绘图参考线

%% 绘图

%% 运动学状态
figure('Name', 'Kinematics', 'Color', 'w');

% figure(1); 
subplot(3, 2, 1); 
plot(t, u, 'LineWidth', 1.5); 
title('纵向速度 u (m/s)'); 
grid on; 
xlabel('时间 (s)');
% ylim([-10, 10]);

% figure(2); 
subplot(3, 2, 2); 
plot(t, v, 'LineWidth', 1.5); 
title('横向速度 v (m/s)'); 
grid on; 
xlabel('时间 (s)');
% ylim([-5, 5]);

% figure(3); 
subplot(3, 2, 3); 
plot(t, abs(1.524 - z_pos), 'LineWidth', 1.5); 
title('垫升高度(m)'); 
grid on; xlabel('时间 (s)');


% figure(4);
subplot(3, 2, 4); 
plot(t, psi_deg, 'LineWidth', 1.5); 
title('艏向角ψ (deg)'); 
grid on; xlabel('时间 (s)');
% ylim([-2, 2]);


% figure(5);
subplot(3, 2, 5); 
plot(t, phi_deg, 'LineWidth', 1.5); 
title('横摇角 Roll (deg)'); grid on; 
xlabel('时间 (s)');
% ylim([-1.5, 1.5]);

% figure(6); 
subplot(3, 2, 6); 
plot(t, theta_deg, 'LineWidth', 1.5); 
title('纵摇角 Pitch (deg)'); 
grid on; 
xlabel('时间 (s)');
% ylim([-2, 2]);

%% 运动轨迹
% figure(7); 
% set(gcf, 'Name', 'Trajectory', 'Color', 'w');
% plot(y_pos, x_pos, 'b-', 'LineWidth', 2); hold on;
% plot(y_pos(1), x_pos(1), 'go', 'MarkerFaceColor', 'g', 'DisplayName', '起点'); 
% plot(y_pos(end), x_pos(end), 'rs', 'MarkerFaceColor', 'r', 'DisplayName', '终点');
% 
% xlabel('东向East(m)'); ylabel('北向North(m)');
% title('气垫船运动轨迹');
% axis equal; grid on; legend show;

figure(7); 
set(gcf, 'Name', 'Trajectory', 'Color', 'w');
plot(y_pos, x_pos, 'b-', 'LineWidth', 2, 'DisplayName', '实际轨迹'); hold on;

% 绘制参考航线 =====
plot(wpt.pos.y, wpt.pos.x, 'k--', 'LineWidth', 1.5, 'DisplayName', '参考路径');
plot(wpt.pos.y, wpt.pos.x, 'ko', 'MarkerFaceColor', 'y', 'DisplayName', '航点');

plot(y_pos(1), x_pos(1), 'go', 'MarkerFaceColor', 'g', 'DisplayName', '起点'); 
plot(y_pos(end), x_pos(end), 'rs', 'MarkerFaceColor', 'r', 'DisplayName', '终点');

xlabel('东向East(m)'); ylabel('北向North(m)');
title('气垫船运动轨迹');
axis equal; grid on; legend show;

%% 侧滑角绘图
% 计算侧滑角 (Sideslip Angle)
beta_rad = atan2(v, u);      % 结果为弧度
beta_deg = rad2deg(beta_rad); % 转换为角度

%  绘制侧滑角图像
figure(8); 
set(gcf, 'Name', 'Sideslip Angle', 'Color', 'w');

plot(t, beta_deg, 'b-', 'LineWidth', 1.5); hold on;
yline(0, 'k--', 'LineWidth', 1, 'DisplayName', '无侧滑');

title('气垫船侧滑角 \beta (Sideslip Angle)');
xlabel('时间 (s)');
ylabel('侧滑角 (deg)');
legend('侧滑角 \beta', 'Location', 'best');
% ylim([-1.5, 1.5]);
grid on;

%% 环境风
figure('Name', 'Wind Environment Verify', 'Color', 'w');

% 存储历史数据
wind_speed_record = zeros(num_steps, 1);
wind_dir_record   = zeros(num_steps, 1); % 新增：用于记录风向

% 获取数据
for k = 1:num_steps
    [wind_speed_record(k), wind_dir_record(k)] = Get_Wind_Step(t(k), wind_data);
end

% 子图1：风速
subplot(2, 1, 1); 
plot(t, wind_speed_record, 'b-', 'LineWidth', 1.2);
title('环境风速');
xlabel('时间(s)'); ylabel('风速(m/s)');
grid on; axis tight;

% 子图2：风向
subplot(2, 1, 2); % 2行1列，第2张图
plot(t, rad2deg(wind_dir_record), 'r-', 'LineWidth', 1.2);
title('环境风向');
xlabel('时间(s)'); ylabel('风向(deg)');
grid on; axis tight;

mean_dir_deg = rad2deg(wind_data.Psi_mean);

%% 气垫压强绘图
figure('Name', 'Cushion Pressure Distribution', 'Color', 'w');

% 转换单位 Pa -> kPa 
P_scale = 1; 

% 前左气室
subplot(2, 2, 1);
plot(t, P_FR * P_scale, 'r-', 'LineWidth', 1.5); hold on;
yline(P_static_si * P_scale, 'k:', 'LineWidth', 1.2);
ylabel('压强 (Pa)'); xlabel('时间 (s)');
title('前左气室 (Front-Right)');
grid on; axis tight;

% 前右气室
subplot(2, 2, 2);
plot(t, P_FL * P_scale, 'b-', 'LineWidth', 1.5); hold on;
yline(P_static_si * P_scale, 'k:', 'LineWidth', 1.2, 'DisplayName', '稳态');
ylabel('压强 (Pa)'); xlabel('时间 (s)');
title('前右气室 (Front-Left)');
grid on; axis tight;


% 后左气室
subplot(2, 2, 3);
plot(t, P_RL * P_scale, 'r-', 'LineWidth', 1.5); hold on;
yline(P_static_si * P_scale, 'k:', 'LineWidth', 1.2);
ylabel('压强 (Pa)'); xlabel('时间 (s)');
title('后左气室 (Rear-Left)');
grid on; axis tight;
% ylim(y_limits);

% 后右气室
subplot(2, 2, 4);
plot(t, P_RR * P_scale, 'b-', 'LineWidth', 1.5); hold on;
yline(P_static_si * P_scale, 'k:', 'LineWidth', 1.2);
ylabel('压强 (Pa)'); xlabel('时间 (s)');
title('后右气室 (Rear-Right)');
grid on; axis tight;
% ylim(y_limits);

% % 后左气室
% subplot(2, 2, 3);
% plot(t, P_RR * P_scale, 'r-', 'LineWidth', 1.5); hold on;
% yline(P_static_si * P_scale, 'k:', 'LineWidth', 1.2);
% ylabel('压强 (Pa)'); xlabel('时间 (s)');
% title('后左气室 (Rear-Left)');
% grid on; axis tight;
% % ylim(y_limits);
% 
% 
% % 后右气室
% subplot(2, 2, 4);
% plot(t, P_RL * P_scale, 'b-', 'LineWidth', 1.5); hold on;
% yline(P_static_si * P_scale, 'k:', 'LineWidth', 1.2);
% ylabel('压强 (Pa)'); xlabel('时间 (s)');
% title('后右气室 (Rear-Right)');
% grid on; axis tight;
% % ylim(y_limits);

grid on;

%% 绘制侧力分量图 ---
% figure('Name', 'Lateral Gravity Force', 'Color', 'w');
% 
% % 需要的参数
% m_slugs = 10879.5;
% SLUG2KG = 14.5939;
% m_kg = m_slugs * SLUG2KG;  % 质量转换
% g_SI = 9.80665;            % 重力加速度
% 
% 
% % phi_rad = sol(:, 4); 
% theta_rad = sol(:, 5); 
% 
% % 计算重力在横向的分量
% % 公式：Fy_g = m * g * sin(phi)
% % Y_G_force = m_kg * g_SI * sin(phi_rad);
% X_G_force = -m_kg * g_SI * sin(theta_rad);
% 
% % 绘图
% % plot(t, Y_G_force, 'b-', 'LineWidth', 1.5);
% plot(t, X_G_force, 'b-', 'LineWidth', 1.5);
% grid on;
% title('重力侧滑分量', 'FontSize', 12);
% xlabel('时间 Time (s)', 'FontSize', 11);
% ylabel('侧向力 Force (N)', 'FontSize', 11);

%% 绘制 MPC 舵角指令图
figure('Name', 'Rudder Command', 'Color', 'w');

plot(t, MPC_rudder_history, 'k-', 'LineWidth', 1.5, 'DisplayName', '舵角指令'); hold on;
% 画出舵角的物理极限线作为参考
yline(30, 'r--', 'LineWidth', 1.2, 'DisplayName', '物理极限');
yline(-30, 'r--', 'LineWidth', 1.2, 'HandleVisibility', 'off');

title('MPC 舵角控制指令 \delta', 'FontSize', 12);
xlabel('时间 Time (s)', 'FontSize', 11);
ylabel('指令舵角 (deg)', 'FontSize', 11);
ylim([-35, 35]); % 稍微留出一点裕度
legend('Location', 'best');
grid on;
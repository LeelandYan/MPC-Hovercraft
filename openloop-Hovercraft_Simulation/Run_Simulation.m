%% 气垫船6DOF仿真 
clear; clc; close all;

%% 初始状态定义
% 状态向量 X = [x, y, z, phi, theta, psi, u, v, w, p, q, r]

x0 = 0;      % 北向位置 (m)
y0 = 0;      % 东向位置 (m)
z0 = 1.5;   % 垂向位移 (m)
phi0 = 0;    % 横摇角 (rad)
theta0 = 0;  % 纵摇角 (rad)
psi0 = deg2rad(0); % 初始艏向 (rad)，0度为正北

% --- 速度 ---
u0 = 0;      % 纵向速度 (m/s)
v0 = 0;      % 横向速度 (m/s)
w0 = 0;      % 垂向速度 (m/s)
p0 = 0;      % 横摇角速度 (rad/s)
q0 = 0;      % 纵摇角速度 (rad/s)
r0 = 0;      % 艏摇角速度 (rad/s)

% 组装初始状态向量
X0 = [x0, y0, z0, phi0, theta0, psi0, u0, v0, w0, p0, q0, r0];

%%
% 设定平均风速，风向
V_wind_mean = 0.0001; 
Psi_wind_mean = deg2rad(0); 
wind_data = Init_Wind(V_wind_mean, Psi_wind_mean);

% 控制指令
rudder_angle = 0; % 舵角指令

%% 使用4阶龙格库塔法
% 定义时间步长和总时长
dt = 0.05;              
T_end = 100;            % 仿真结束时间
t = 0:dt:T_end;         % 生成时间向量
num_steps = length(t);  % 总步数

% 初始化存储数组
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
    
    % --- RK4迭代 ---
    [dX1, ~] = model_jeff_b(t_curr,          X_curr,              rudder_angle, wind_data);
    [dX2, ~] = model_jeff_b(t_curr + 0.5*dt, X_curr + 0.5*dt*dX1, rudder_angle, wind_data);
    [dX3, ~] = model_jeff_b(t_curr + 0.5*dt, X_curr + 0.5*dt*dX2, rudder_angle, wind_data);
    [dX4, ~] = model_jeff_b(t_curr + dt,     X_curr + dt*dX3,     rudder_angle, wind_data);
    X_next = X_curr + (dt / 6) * (dX1 + 2*dX2 + 2*dX3 + dX4);
   
    sol(k+1, :) = X_next';
    
    [~, P_out] = model_jeff_b(t(k+1), X_next, rudder_angle, wind_data);
    P_hist_Pa(k+1, :) = P_out(:)';
    
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
figure(7); 
set(gcf, 'Name', 'Trajectory', 'Color', 'w');
plot(y_pos, x_pos, 'b-', 'LineWidth', 2); hold on;
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
figure('Name', 'Lateral Gravity Force', 'Color', 'w');

% 需要的参数
m_slugs = 10879.5;
SLUG2KG = 14.5939;
m_kg = m_slugs * SLUG2KG;  % 质量转换
g_SI = 9.80665;            % 重力加速度

% 提取数据,sol(:,4)是横摇角
phi_rad = sol(:, 4); 

% 计算重力在横向的分量
% 公式：Fy_g = m * g * sin(phi)
Y_G_force = m_kg * g_SI * sin(phi_rad);

% 绘图
plot(t, Y_G_force, 'b-', 'LineWidth', 1.5);
grid on;
title('重力侧滑分量', 'FontSize', 12);
xlabel('时间 Time (s)', 'FontSize', 11);
ylabel('侧向力 Force (N)', 'FontSize', 11);

% -------------------------------------------------------------------------








% -------------------------------------------------------------------------
% %% 打舵回转对比
% clear; clc; close all;
% 
% %% 仿真设置
% % 定义要对比的舵角列表
% rudder_test_values = [5, 10, 15]; 
% line_colors = {'b', 'b', 'b'};     % 对应颜色：全为蓝色
% % --- 新增：定义线型列表 ---
% % '-' = 实线, '--' = 虚线, ':' = 点线
% line_styles = {'-', '--', ':'};    
% 
% % 初始化风 ---
% V_wind_mean = 3;              % 平均风速 (m/s)
% Psi_wind_mean = deg2rad(90);   % 平均风向 (rad)
% % 请确保您的目录下有 Init_Wind 函数
% wind_data = Init_Wind(V_wind_mean, Psi_wind_mean); 
% 
% % 定义时间步长和总时长
% dt = 0.1;               
% T_end = 400;            
% t = 0:dt:T_end;          
% num_steps = length(t);  
% 
% %% 准备绘图
% % 轨迹对比
% fig_traj = figure('Name', 'Trajectory Comparison', 'Color', 'w');
% hold on; grid on; axis equal;
% xlabel('东向 East (m)'); ylabel('北向 North (m)');
% title('气垫船运动轨迹');
% 
% % 侧滑角对比
% fig_beta = figure('Name', 'Sideslip Angle Comparison', 'Color', 'w');
% hold on; grid on;
% xlabel('时间 (s)'); ylabel('侧滑角 (deg)');
% title('侧滑角');
% 
% % 纵向速度对比
% fig_u = figure('Name', 'Surge Velocity Comparison', 'Color', 'w');
% hold on; grid on;
% xlabel('时间 (s)'); ylabel('纵向速度 u (m/s)');
% title('纵向速度');
% 
% % 横向速度对比
% fig_v = figure('Name', 'Sway Velocity Comparison', 'Color', 'w');
% hold on; grid on;
% xlabel('时间 (s)'); ylabel('横向速度 v (m/s)');
% title('横向速度');
% 
% % 艏向角对比图
% fig_psi = figure('Name', 'Heading Angle Comparison', 'Color', 'w');
% hold on; grid on;
% xlabel('时间 (s)'); ylabel('艏向角 \psi (deg) [0-360]');
% title('艏向角对比 (Heading)');
% ylim([0, 360]); % 锁定Y轴范围，方便观察罗盘方位
% yticks(0:45:360); % 设置刻度间隔，符合航海习惯
% 
% %% 循环
% for i = 1:length(rudder_test_values)
%     
%     % --- 获取当前工况参数 ---
%     current_angle = rudder_test_values(i);
%     current_color = line_colors{i};
%     % --- 新增：获取当前线型 ---
%     current_style = line_styles{i};
% 
%     fprintf('正在进行舵角 = %d 度的仿真...\n', current_angle);
%     
%     clear model_jeff_b; 
%     
%     % --- 初始状态定义 (每次都要重置) ---
%     x0 = 0; y0 = 0; z0 = 1.5;   
%     phi0 = 0; theta0 = 0; psi0 = deg2rad(0); 
%     u0 = 0; v0 = 0; w0 = 0;     
%     p0 = 0; q0 = 0; r0 = 0;     
%     X0 = [x0, y0, z0, phi0, theta0, psi0, u0, v0, w0, p0, q0, r0];
%     
%     % 初始化存储数组
%     sol = zeros(num_steps, 12); 
%     sol(1, :) = X0;          
% 
%     % --- 龙格库塔循环 ---
%     for k = 1 : num_steps - 1
%             t_curr = t(k);
%             X_curr = sol(k, :)'; 
%             
%             % 传递 wind_data 参数
%             [dX1, ~] = model_jeff_b(t_curr,          X_curr,               current_angle, wind_data);
%             [dX2, ~] = model_jeff_b(t_curr + 0.5*dt, X_curr + 0.5*dt*dX1, current_angle, wind_data);
%             [dX3, ~] = model_jeff_b(t_curr + 0.5*dt, X_curr + 0.5*dt*dX2, current_angle, wind_data);
%             [dX4, ~] = model_jeff_b(t_curr + dt,     X_curr + dt*dX3,     current_angle, wind_data);
%             
%             X_next = X_curr + (dt / 6) * (dX1 + 2*dX2 + 2*dX3 + dX4);
%             sol(k+1, :) = X_next';
%      end
%     
%     %% 数据解包
%     x_pos = sol(:,1);
%     y_pos = sol(:,2);
%     u_vel = sol(:,7);
%     v_vel = sol(:,8);
% 
%     % 艏向角处理 
%     psi_rad = sol(:,6); 
%     psi_deg_raw = rad2deg(psi_rad);      % 原始角度
%     psi_deg_360 = mod(psi_deg_raw, 360); % 取模运算，映射到 0-360
%     
%     % 计算侧滑角
%     beta_deg = rad2deg(atan2(v_vel, u_vel));
% 
%     % --- 绘图部分 (增加了 LineStyle 和调整了 LineWidth) ---
%     
%     % 绘制轨迹 
%     figure(fig_traj);
%     plot(y_pos, x_pos, 'LineWidth', 1.5, 'LineStyle', current_style, 'Color', current_color, ...
%          'DisplayName', sprintf('舵角 %d^{\\circ}', current_angle));
%     % 标记终点
%     plot(y_pos(end), x_pos(end), 's', 'MarkerFaceColor', current_color, 'HandleVisibility', 'off');
% 
%     % 绘制侧滑角
%     figure(fig_beta);
%     plot(t, beta_deg, 'LineWidth', 1.5, 'LineStyle', current_style, 'Color', current_color, ...
%          'DisplayName', sprintf('舵角 %d^{\\circ}', current_angle));
% 
%     % 纵向速度
%     figure(fig_u);
%     plot(t, u_vel, 'LineWidth', 1.5, 'LineStyle', current_style, 'Color', current_color, ...
%          'DisplayName', sprintf('舵角 %d^{\\circ}', current_angle));
% 
%     % 横向速度
%     figure(fig_v);
%     plot(t, v_vel, 'LineWidth', 1.5, 'LineStyle', current_style, 'Color', current_color, ...
%          'DisplayName', sprintf('舵角 %d^{\\circ}', current_angle));
% 
%     % 艏向角
%     figure(fig_psi);
%     plot(t, psi_deg_360, 'LineWidth', 1.5, 'LineStyle', current_style, 'Color', current_color, ...
%          'DisplayName', sprintf('舵角 %d^{\\circ}', current_angle));
% 
% end
% 
% % 轨迹图装饰
% figure(fig_traj);
% plot(0, 0, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', '起点', 'MarkerSize', 10); % 起点加大一点
% legend('show', 'Location', 'best');    % 显示图例
% 
% % 其他图装饰
% figure(fig_beta);
% yline(0, 'k--', 'HandleVisibility', 'off');
% legend('show', 'Location', 'best');
% 
% figure(fig_u);
% legend('show', 'Location', 'best');
% 
% figure(fig_v);
% yline(0, 'k--', 'HandleVisibility', 'off');
% legend('show', 'Location', 'best');
% 
% figure(fig_psi);
% legend('show', 'Location', 'best');
% 
% fprintf('仿真完成。\n');

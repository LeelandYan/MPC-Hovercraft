function [V_now, Psi_now] = Get_Wind_Step(time, wind_param)
    % Get_Wind_Step 计算特定时间点的风速和风向 
    
    % 提取参数
    omega = wind_param.omega;
    V_mean = wind_param.V_mean;
    Psi_mean = wind_param.Psi_mean;
    
    % 1. 计算纵向湍流速度 u_gust
    u_gust = sum(wind_param.amp_long .* cos(omega * time + wind_param.epsilon_long));
    
    % 2. 计算侧向湍流速度 v_gust 
    v_gust = sum(wind_param.amp_lat  .* cos(omega * time + wind_param.epsilon_lat));
    
    % 3. 合成
    V_total_x = V_mean + u_gust;
    V_total_y = v_gust;
    
    % 4. 
    % 瞬时风速大小
    V_now = sqrt(V_total_x^2 + V_total_y^2);
    
    % 瞬时风向角度 = 平均风向 + 扰动角度
    angle_perturbation = atan2(V_total_y, V_total_x);
    Psi_now = Psi_mean + angle_perturbation;
    
end
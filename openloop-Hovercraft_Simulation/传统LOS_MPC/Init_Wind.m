function wind_param = Init_Wind(V_mean, Psi_mean)
    % Init_Wind 初始化风场参数 
    % 输入：平均风速 V_mean, 平均风向 Psi_mean
    % 输出：结构体 wind_param
    
    % 1. 基础参数
    z = 10; % 高度 10m
    if(z <= 20)
        alpha = -0.125;
    else
        alpha = -0.275;
    end
    
    % 2. 谱参数计算 (基于绝对风速)
    % 纵向湍流强度
    sigma_u = 0.15 * (z/20)^alpha * V_mean; 
    omega_p = 2 * pi * 0.0025 * V_mean;
    
    % 3. 生成频率范围
    N = 50; 
    omega = linspace(0.001, 2, N); 
    delta_omega = omega(2) - omega(1);
    
    % 4. 计算功率谱 
    % S_long: 纵向谱 
    S_long = sigma_u^2 ./ (omega_p * (1 + 1.5*omega./omega_p).^(5.0/3.0));
    
    % S_lat: 侧向谱
    S_lat = 0.1 * S_long; 
    
    % 5. 计算振幅
    amp_long = sqrt(2 * S_long * delta_omega);
    amp_lat  = sqrt(2 * S_lat  * delta_omega);
    
    % 6. 生成两组独立的随机相位 
    % epsilon_long 用于风速波动
    % epsilon_lat  用于风向波动
    epsilon_long = rand(1, N) * 2 * pi; 
    epsilon_lat  = rand(1, N) * 2 * pi; 
    
    % 7. 打包数据
    wind_param.V_mean = V_mean;
    wind_param.Psi_mean = Psi_mean;
    wind_param.omega = omega;
    wind_param.amp_long = amp_long;     % 纵向振幅
    wind_param.amp_lat  = amp_lat;      % 侧向振幅
    wind_param.epsilon_long = epsilon_long;
    wind_param.epsilon_lat  = epsilon_lat;
    wind_param.N = N;
end
function [psi_ref, y_e, k_active] = LOSpsi(x, y, Delta_h, R_switch, wpt)
% 基础LOS
%
% 输入:   
%   x, y:     当前的北东坐标 (m)
%   Delta_h:  前视距离常数 (m)
%   R_switch: 航路点切换半径 (m)
%   wpt:      结构体，包含航路点坐标 wpt.pos.x 和 wpt.pos.y
%
% 输出:  
%   psi_ref:  期望航向角 (rad)
%   y_e:      横向偏迹误差 (m)
%   k_active: 当前正在追踪的目标航路点索引

persistent k;        % 活跃航路点索引
persistent xk yk;    % 活跃航路点坐标


if isempty(k)
    % 检查 R_switch 是否合理 (必须小于航路点间最小距离)
    if R_switch > min( sqrt( diff(wpt.pos.x).^2 + diff(wpt.pos.y).^2 ) )
        error("航路点之间的距离必须大于切换半径 R_switch");
    end
    if (R_switch < 0); error("R_switch 必须大于 0"); end
    if (Delta_h < 0); error("Delta_h 必须大于 0"); end

    k = 1;                  % 初始化索引为第1个航路点
    xk = wpt.pos.x(k);
    yk = wpt.pos.y(k);
end


% 读取下一个航路点 (xk_next, yk_next)
n = length(wpt.pos.x);
if k < n                        % 如果还有下一个航路点，读取它
    xk_next = wpt.pos.x(k+1);  
    yk_next = wpt.pos.y(k+1);   
else                            % 如果已经是最后一个航路点，保持最后一段航向无限延伸
    bearing = atan2((wpt.pos.y(n)- wpt.pos.y(n-1)), (wpt.pos.x(n)-wpt.pos.x(n-1)));
    R = 1e10;
    xk_next = wpt.pos.x(n) + R * cos(bearing);
    yk_next = wpt.pos.y(n) + R * sin(bearing); 
end


% 计算路径参数与误差
% 计算当前路径段的切向角 pi_h (相对于北向)
pi_h = atan2( (yk_next-yk), (xk_next-xk) ); 

% 计算纵向误差 x_e 和横向误差 y_e (坐标系旋转)
x_e =  (x-xk) * cos(pi_h) + (y-yk) * sin(pi_h);
y_e = -(x-xk) * sin(pi_h) + (y-yk) * cos(pi_h);

% 航路点切换
d = sqrt( (xk_next-xk)^2 + (yk_next-yk)^2 );
if ( (d - x_e < R_switch) && (k < n) )
    k = k + 1;
    xk = xk_next;       % 更新当前活跃航路点
    yk = yk_next; 
end


% LOS制导
psi_ref = pi_h - atan( y_e/Delta_h ); 

k_active = k; % 返回当前活跃的航路点索引

end
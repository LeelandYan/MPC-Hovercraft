function [psi_ref,y_e,beta_c_hat,k_active] = ALOSpsi(x,y,Delta_h,gamma_h,h,R_switch,wpt)

persistent k;        % active waypoint index, initialized by 'clear ALOSpsi'
persistent xk yk;    % active waypoint (xk, yk) corresponding to integer k
persistent beta_hat; % estimate of the crab angle

% ------------------------------------------------------------------------------
% Initialization of (xk, yk) and (xk_next, yk_next), and integral state 
% ------------------------------------------------------------------------------
if isempty(k)

    % Check if R_switch is smaller than the minimum distance between the waypoints
    if R_switch > min( sqrt( diff(wpt.pos.x).^2 + diff(wpt.pos.y).^2 ) )
        error("The distances between the waypoints must be larger than R_switch");
    end

    % Check input parameters
    if (R_switch < 0); error("R_switch must be larger than zero"); end
    if (Delta_h < 0); error("Delta must be larger than zero"); end

    beta_hat = 0;           % initial states
    k = 1;                  % set first waypoint as the active waypoint
    xk = wpt.pos.x(k);
    yk = wpt.pos.y(k);
    fprintf('Active waypoints:\n')
    fprintf('  (x%1.0f, y%1.0f) = (%.2f, %.2f) \n',k,k,xk,yk);

end

% ------------------------------------------------------------------------------
% Read next waypoint (xk_next, yk_next) from wpt.pos 
% ------------------------------------------------------------------------------
n = length(wpt.pos.x);
if k < n                        % if there are more waypoints, read next one 
    xk_next = wpt.pos.x(k+1);  
    yk_next = wpt.pos.y(k+1);    
else                            % else, continue with last bearing
    bearing = atan2((wpt.pos.y(n)- wpt.pos.y(n-1)), (wpt.pos.x(n)-wpt.pos.x(n-1)));
    R = 1e10;
    xk_next = wpt.pos.x(n) + R * cos(bearing);
    yk_next = wpt.pos.y(n) + R * sin(bearing); 
end

% ------------------------------------------------------------------------------
% Compute the path-tangential angle w.r.t. north
% ------------------------------------------------------------------------------
pi_h = atan2( (yk_next-yk), (xk_next-xk) ); 

% Along-track and cross-track errors (x_e, y_e) 
x_e =  (x-xk) * cos(pi_h) + (y-yk) * sin(pi_h);
y_e = -(x-xk) * sin(pi_h) + (y-yk) * cos(pi_h);

% If the next waypoint satisfy the switching criterion, k = k + 1
d = sqrt( (xk_next-xk)^2 + (yk_next-yk)^2 );
if ( (d - x_e < R_switch) && (k < n) )
    k = k + 1;
    xk = xk_next;       % update active waypoint
    yk = yk_next; 
    fprintf('  (x%1.0f, y%1.0f) = (%.2f, %.2f) \n',k,k,xk,yk);
end

% ------------------------------------------------------------------------------
% ALOS guidance law
% ------------------------------------------------------------------------------
psi_ref = pi_h - beta_hat - atan( y_e/Delta_h ); 

beta_c_hat = beta_hat; % Return crab angle estimate
k_active = k; % Return active waypoint index

% Propagation of crab angle estimate to time k+1
beta_hat = beta_hat + h * gamma_h * Delta_h * y_e / sqrt( Delta_h^2 + y_e^2 );

% ===== 新增：抗积分饱和 (Anti-Windup) 限幅 =====
beta_max = deg2rad(20); % 限制最大侧滑角估计值为 20 度
if beta_hat > beta_max
    beta_hat = beta_max;
elseif beta_hat < -beta_max
    beta_hat = -beta_max;
end

end

% -------------------------------------------------------------------------
% MATLAB Code: Digital Twin Control - Case 2 (Perfect Initialization)
% -------------------------------------------------------------------------

clear; clc; close all;

%% 1. 初始化与参数
N_sim = 200;            
N_pred = 10;            

% --- 控制器参数 (Eq 63) ---
hp = 1.0; Kpp = 0.75; Kpi = 0.035;            

% --- 模型参数 ---
theta_p = [0.5; 1.5; 0.3];                     % 真实非线性对象参数
theta_d = [1.3178; -0.3866; 1.4238; -1.2606];  % 数字孪生线性参数

% --- 数据存储 ---
% 为了方便对应原图的 t=0..200，我们增加一点缓冲
padding = 20; 
total_len = N_sim + padding + N_pred; 

% 真实系统状态
y_p = zeros(total_len, 1);   
u_p = zeros(total_len, 1);   
z_p = zeros(total_len, 1);   
r   = zeros(total_len, 1);   

% 预测数据存储 (初始化为NaN，方便看出哪里没数据)
y_d_pred_10step = nan(total_len, 1); 
u_d_pred_10step = nan(total_len, 1); 

% --- 生成参考信号 r(t) ---
% t=0对应数组索引 padding+1
for k = 1:total_len
    t_sim = k - padding - 1; % 使得 t_sim 从 0 开始
    if (t_sim >= 0 && t_sim < 50) || (t_sim >= 100 && t_sim < 150)
        r(k) = 1.0;
    else
        r(k) = 0.0;
    end
end

%% 2. 主仿真循环
start_idx = padding + 1;       % 对应 t=0
end_idx   = padding + N_sim;   % 对应 t=200

for k = start_idx : end_idx
    
    % 当前时刻的索引 k
    % 注意：MATLAB索引从1开始，物理时间 t = k - start_idx
    
    % --- A. 实时控制 (Real-time Control) ---
    % 参数提取
    thd1 = theta_d(1); thd2 = theta_d(2); thd3 = theta_d(3); thd4 = theta_d(4);
    
    % 1. 计算控制律 u_p(k)
    % 使用上一时刻的状态 k-1 (对应物理时刻 t-1)
    % 注意：Eq 61 中的 y(t) 在代码循环中其实是已知的当前测量值 y_p(k) (如果在更新方程之后)
    % 但在标准离散步中，通常顺序是：测量y(k) -> 计算u(k) -> 系统演变y(k+1)
    
    % 计算 u(t)
    term_y      = (Kpp - thd1) * y_p(k);
    term_y_prev = -thd2 * y_p(k-1);
    term_z      = Kpi * z_p(k);
    term_u_prev = -thd4 * u_p(k-1);
    
    u_p(k) = (term_y + term_y_prev + term_z + term_u_prev) / thd3;
    
    % 2. 积分更新 z(t+1)
    z_p(k+1) = z_p(k) + hp * (r(k) - y_p(k));
    
    % --- B. 真实对象演变 (Nonlinear Plant) ---
    % y(t+1) calculation
    num = theta_p(1) * y_p(k) * y_p(k-1);
    den = 1 + y_p(k-1)^2;
    nonlinear_y = num / den;
    control_u   = theta_p(2) * u_p(k) + theta_p(3) * u_p(k-1);
    
    y_p(k+1) = nonlinear_y + control_u;
           
    % --- C. 10步超前预测 (Prediction Logic) ---
    % 我们站在时刻 k，要把 k+1 到 k+10 的事情都预测出来
    
    % 初始化虚拟状态 (从当前真实状态出发)
    y_curr = y_p(k);      
    y_prev = y_p(k-1);
    u_curr = u_p(k);
    u_prev = u_p(k-1);
    z_curr = z_p(k+1); % 下一步的积分状态起点
    
    % 临时轨迹容器
    traj_y = zeros(N_pred, 1);
    traj_u = zeros(N_pred, 1);
    
    for i = 1:N_pred
        idx_future = k + i; % 未来的绝对索引
        r_future = r(idx_future);
        
        % 1. 预测模型 (Linear Digital Twin)
        y_next = thd1*y_curr + thd2*y_prev + thd3*u_curr + thd4*u_prev;
        
        % 2. 虚拟控制律
        u_term_y      = (Kpp - thd1) * y_next;
        u_term_y_prev = -thd2 * y_curr;
        u_term_z      = Kpi * z_curr;
        u_term_u_prev = -thd4 * u_curr;
        
        u_next = (u_term_y + u_term_y_prev + u_term_z + u_term_u_prev) / thd3;
        
        % 3. 虚拟积分更新
        z_next = z_curr + hp * (r_future - y_next);
        
        % 状态滚动
        y_prev = y_curr; y_curr = y_next;
        u_prev = u_curr; u_curr = u_next;
        z_curr = z_next;
        
        traj_y(i) = y_next;
        traj_u(i) = u_next;
    end
    
    % --- 关键修改：数据填充策略 ---
    
    % 策略1：标准存储 (保存 t+10 的预测值)
    if k + N_pred <= total_len
        y_d_pred_10step(k + N_pred) = traj_y(N_pred);
        u_d_pred_10step(k + N_pred) = traj_u(N_pred);
    end
    
    % 策略2：初始化填充 (解决前10步空白问题)
    % 如果这是仿真的第一步 (t=0)，那么这一步预测出来的未来10步轨迹，
    % 就是系统在前10步唯一的"预测参考"。
    if k == start_idx
        for i = 1:N_pred
            y_d_pred_10step(k+i) = traj_y(i);
            u_d_pred_10step(k+i) = traj_u(i);
        end
        % 甚至把当前时刻也填上，保证从0开始连线
        y_d_pred_10step(k) = y_p(k);
        u_d_pred_10step(k) = u_p(k);
    end
end

%% 3. 绘图 (Figure 9 & 10)

% 提取绘图范围 (0 到 200)
plot_range = start_idx : end_idx;
t_axis = 0 : (length(plot_range)-1); % 时间轴 0,1,2...200

r_data = r(plot_range);
y_data = y_p(plot_range);
u_data = u_p(plot_range);
ypred_data = y_d_pred_10step(plot_range);
upred_data = u_d_pred_10step(plot_range);

% --- Figure 9: Outputs ---
figure(1); set(gcf, 'Color', 'w', 'Position', [100 100 600 400]);
plot(t_axis, r_data, 'k--', 'LineWidth', 1.2); hold on;
plot(t_axis, y_data, 'b-', 'LineWidth', 1.5);
% 红���预测线
plot(t_axis, ypred_data, 'r-', 'LineWidth', 1.5);

title('Fig. 9 Outputs of the system (Case 2)');
xlabel('Time (steps)'); ylabel('Outputs');
legend('Reference', 'Output y_p(t)', 'Prediction y_d(t|t-10)', 'Location', 'SouthEast');
grid on; axis([0 200 -0.2 1.4]);

% --- Figure 10: Control Inputs ---
figure(2); set(gcf, 'Color', 'w', 'Position', [750 100 600 400]);
% 蓝色控制线
plot(t_axis, u_data, 'b-', 'LineWidth', 1.5); hold on;
% 红色预测线
plot(t_axis, upred_data, 'r-', 'LineWidth', 1.5);

title('Fig. 10 Control inputs of the system (Case 2)');
xlabel('Time (steps)'); ylabel('Control Inputs');
legend('Control u_p(t)', 'Prediction u_d(t|t-10)', 'Location', 'SouthEast');
grid on; axis([0 200 -0.5 0.5]); % 调整纵坐标范围以匹配原图细节

% 为了看清前10步，可以在这里暂停或放大
% xlim([0 50]); 
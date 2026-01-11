%% 数字孪生控制系统 - 案例1: 自适应模型跟踪 (最终修正版)
% 修复内容：
% 1. 循环从 k=1 (t=0) 开始，捕捉初始控制信号
% 2. 修正 PI 控制器时序 (先输出后积分)
% 3. 配色方案与文献 Figure 8 完全一致
% 环境：MATLAB 2025a

clear; clc; close all;

%% 1. 仿真参数
T_sim = 200;
t_vec = 1:T_sim; % 对应 t=0 到 t=199

% 参考输入 r(t)
r = zeros(1, T_sim);
for k = 1:T_sim
    if (k <= 50) || (k > 100 && k <= 150)
        r(k) = 1;
    else
        r(k) = 0;
    end
end

%% 2. 初始化
yp = zeros(1, T_sim + 1); % 多给一个空间防止溢出
up = zeros(1, T_sim + 1);
zp = zeros(1, T_sim + 1); 

yd = zeros(1, T_sim + 1);
ud = zeros(1, T_sim + 1);
zd = zeros(1, T_sim + 1);

% PI 参数 (Eq 63)
K_p = 0.75;
K_i = 0.035;

% RLS 初始化
n_params = 4;
% theta_hat = [th1; th2; th3; th4]
theta_hat = [1.3178; -0.3866; 1.4238; -1.2606]; 
P_cov = 10^8 * eye(n_params); % 初始大协方差
lambda = 0.95; 

% 历史记录
theta_history = zeros(n_params, T_sim);
theta_history(:,1) = theta_hat; % 记录初始值

%% 3. 全过程仿真循环 (从 t=0 开始)
fprintf('开始高精度复现仿真...\n');

for k = 1:T_sim
    % 当前时刻 t = k-1
    % 我们需要利用 t-1, t-2 的数据计算 t 时刻的控制和输出
    
    % --- Part A: 设定真实对象参数 ---
    if k <= 100
        theta_p1 = 0.5; theta_p2 = 1.5; theta_p3 = 0.3;
    else
        theta_p1 = 0.8; theta_p2 = 0.7; theta_p3 = 0.5;
    end
    
    % --- Part B: RLS 自适应 (在计算控制律之前更新模型) ---
    % 构建回归向量 phi(t-1) = [y(t-1); y(t-2); u(t-1); u(t-2)]
    % 处理边界条件：如果下标 < 1，则值为 0
    y_tm1 = get_val(yp, k-1); y_tm2 = get_val(yp, k-2);
    u_tm1 = get_val(up, k-1); u_tm2 = get_val(up, k-2);
    
    phi = [y_tm1; y_tm2; u_tm1; u_tm2];
    
    % 只有当 k > 1 (即 t > 0) 且有数据时才真正更新，但协方差矩阵 P 会随遗忘因子增长
    % 真实测量值 y_p(t) 即 yp(k)
    % 注意：yp(k) 是在上一循环步 (k-1) 的末尾计算出来的，代表当前时刻 t 的状态
    y_meas = yp(k); 
    
    % RLS Update
    % 计算增益
    K_gain = (P_cov * phi) / (lambda + phi' * P_cov * phi);
    
    % 计算预测误差 (Prior Error)
    error_est = y_meas - phi' * theta_hat;
    
    % 更新参数
    theta_hat = theta_hat + K_gain * error_est;
    
    % 更新协方差
    P_cov = (eye(n_params) - K_gain * phi') * P_cov / lambda;
    
    % 记录
    theta_history(:, k) = theta_hat;
    
    % --- Part C: 数字控制器 D-Controller ---
    % 1. 计算控制律 u_d(t) (Eq 59)
    % 使用当前积分状态 zd(k) (即 z_d(t))
    % 以及 t-1 时刻的 y 和 u
    y_dtm1 = get_val(yd, k-1);
    u_dtm1 = get_val(ud, k-1);
    
    th1 = theta_hat(1); th2 = theta_hat(2);
    th3 = theta_hat(3); th4 = theta_hat(4);
    
    if abs(th3) < 1e-4, th3 = 1e-4 * sign(th3); end % 奇异值保护
    
    % u_d(t)
    num_d = (K_p - th1)*yd(k) - th2*y_dtm1 + K_i*zd(k) - th4*u_dtm1;
    ud(k) = num_d / th3;
    
    % 2. 更新积分状态 z_d(t+1) (Eq 58)
    zd(k+1) = zd(k) + (r(k) - yd(k));
    
    % 3. 数字模型预测下一时刻输出 y_d(t+1)
    % phi_next = [y(t); y(t-1); u(t); u(t-1)]
    phi_next = [yd(k); y_dtm1; ud(k); u_dtm1];
    yd(k+1) = theta_hat' * phi_next;
    
    % --- Part D: 真实控制器 R-Controller & Plant ---
    % 1. 计算控制律 u_p(t) (Eq 61)
    % 参数与 D-Controller 相同
    num_p = (K_p - th1)*yp(k) - th2*y_tm1 + K_i*zp(k) - th4*u_tm1;
    up(k) = num_p / th3;
    
    % 2. 更新积分状态 z_p(t+1) (Eq 60)
    zp(k+1) = zp(k) + (r(k) - yp(k));
    
    % 3. 真实对象非线性动力学 y_p(t+1)
    % y_p(t+1) = (theta_p1 * y(t) * y(t-1)) / (1 + y(t-1)^2) + ...
    term_nl = (theta_p1 * yp(k) * y_tm1) / (1 + y_tm1^2);
    yp(k+1) = term_nl + theta_p2 * up(k) + theta_p3 * u_tm1;
    
end

%% 4. 绘图 (匹配文献配色)
figure('Color','w', 'Position', [100, 100, 800, 800]);

% 子图1: 输出响应
subplot(2,1,1);
plot(t_vec, r(1:T_sim), 'k--', 'LineWidth', 1.5); hold on;
plot(t_vec, yp(1:T_sim), 'b-', 'LineWidth', 1.5); 
plot(t_vec, yd(1:T_sim), 'r:', 'LineWidth', 2);
line([100 100], [-0.5 1.5], 'Color', 'b', 'LineStyle', '-');
text(102, 0.5, '\leftarrow Plant Change', 'Color', 'b', 'FontSize', 10);
title('Output Response');
legend('Reference', 'Real Plant', 'Digital Twin');
grid on; ylim([-0.2, 1.2]);

% 子图2: 参数估计 (Fig 8 复现)
subplot(2,1,2);
% 配色方案 (参考 Fig 8):
% Theta1: Blue (蓝色)
% Theta2: Red/Orange (红色/红褐色)
% Theta3: Yellow/Orange (黄色)
% Theta4: Purple (紫色)
color_th1 = [0, 0.4470, 0.7410];      % 蓝
color_th2 = [0.8500, 0.3250, 0.0980]; % 红/橙
color_th3 = [0.9290, 0.6940, 0.1250]; % 黄
color_th4 = [0.4940, 0.1840, 0.5560]; % 紫

plot(t_vec, theta_history(1,:), 'Color', color_th1, 'LineWidth', 1.5); hold on;
plot(t_vec, theta_history(2,:), 'Color', color_th2, 'LineWidth', 1.5);
plot(t_vec, theta_history(3,:), 'Color', color_th3, 'LineWidth', 1.5);
plot(t_vec, theta_history(4,:), 'Color', color_th4, 'LineWidth', 1.5);

% 标记参数突变点
line([100 100], [-4 4], 'Color', 'b', 'LineStyle', '-', 'LineWidth', 2);
text(70, -3.5, 'The plant parameters change \uparrow', 'Color', 'k', 'FontSize', 10, 'FontWeight', 'bold');

title('Parameters of the digital plant (Fig 8 Replication)');
legend('\theta_{d1}(t)', '\theta_{d2}(t)', '\theta_{d3}(t)', '\theta_{d4}(t)', 'Location', 'Best');
xlabel('Time (steps)'); ylabel('Parameter Value');
grid on; 
ylim([-4, 3]); % 匹配文献纵坐标范围

%% 辅助函数: 安全获取数组值
function val = get_val(arr, idx)
    if idx < 1
        val = 0;
    else
        val = arr(idx);
    end
end
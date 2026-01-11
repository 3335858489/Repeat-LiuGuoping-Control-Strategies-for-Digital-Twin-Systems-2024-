% 清除工作区
clear; clc; close all;

%% 1. 初始化设置
T_end = 200;                    % 仿真时间 0 到 200
time = 0:T_end;                 % 时间向量，长度 201
len_plot = length(time);        % 绘图点数

% 物理对象参数
theta_p = [0.5; 1.5; 0.3];

% 数字孪生参数 (初始估计值)
theta_d = [1.3178; -0.3866; 1.4238; -1.2606];

% 控制器参数
% 初始参数 (t < 100)
K_p_init = [0.75, 0.035];       % [Kpp/Kdp, Kpi/Kdi]
% 优化参数 (t >= 100)
K_p_opt  = [0.7, 0.04];         % [Kpp/Kdp, Kpi/Kdi]

%% 2. 数据数组预分配 (含历史数据缓冲)
% 为了处理 t-1, t-2，我们在数组前面加 2 个缓冲位
offset = 2;
total_len = len_plot + offset; 

% 初始化所有信号为 0
r_all  = zeros(1, total_len);   % 方波参考信号
yr_all = zeros(1, total_len);   % 期望平滑输出 (Eq. 62)
yp_all = zeros(1, total_len);   % 物理输出
up_all = zeros(1, total_len);   % 物理控制输入
zp_all = zeros(1, total_len);   % 物理积分状态

yd_all = zeros(1, total_len);   % 数字孪生输出
ud_all = zeros(1, total_len);   % 数字孪生控制输入
zd_all = zeros(1, total_len);   % 数字孪生积分状态

%% 3. 仿真主循环
% 循环变量 t 代表物理时间 0, 1, ..., 200
for t = 0 : T_end
    
    % k 是数组中的索引 (因为有 offset)
    k = t + 1 + offset; 
    
    % --- (1) 生成方波信号 r(t) ---
    if (t >= 0 && t < 50) || (t >= 100 && t < 150)
        r_val = 1;
    else
        r_val = 0;
    end
    r_all(k) = r_val;
    
    % --- (2) 计算期望输出 yr(t) (Eq. 62) ---
    % 注意：yr 仅作为参考曲线绘制，不参与控制器的反馈计算
    % yr(t) = 1.7210*yr(t-1) - 0.7558*yr(t-2) + 0.0182*r(t-1) + 0.0166*r(t-2)
    yr_all(k) = 1.7210 * yr_all(k-1) - 0.7558 * yr_all(k-2) + ...
                0.0182 * r_all(k-1) + 0.0166 * r_all(k-2);

    % --- (3) 参数切换 (Case 3 核心) ---
    if t < 100
        Kp = K_p_init(1); Ki = K_p_init(2);
    else
        Kp = K_p_opt(1);  Ki = K_p_opt(2);
    end
    
    % --- (4) 控制器计算 ---
    % 注意：积分状态 zp(k), zd(k) 是基于上一步更新得来的，代表当前时刻 t 的积分值
    
    % >> 数字孪生控制器 D-controller (Eq. 59)
    % ud(t) = ...
    num_d = (Kp - theta_d(1))*yd_all(k) - theta_d(2)*yd_all(k-1) + Ki*zd_all(k) - theta_d(4)*ud_all(k-1);
    ud_all(k) = num_d / theta_d(3);
    
    % >> 物理控制器 R-controller (Eq. 61)
    % up(t) = ...
    num_p = (Kp - theta_d(1))*yp_all(k) - theta_d(2)*yp_all(k-1) + Ki*zp_all(k) - theta_d(4)*up_all(k-1);
    up_all(k) = num_p / theta_d(3);
    
    % 饱和限制 (Optional)
    if up_all(k) > 10, up_all(k)=10; elseif up_all(k)<-10, up_all(k)=-10; end
    if ud_all(k) > 10, ud_all(k)=10; elseif ud_all(k)<-10, ud_all(k)=-10; end
    
    % --- (5) 系统演化 (计算下一时刻 t+1 的输出) ---
    
    % 物理对象非线性方程
    % yp(t+1) = ...
    term1 = (theta_p(1) * yp_all(k) * yp_all(k-1)) / (1 + yp_all(k-1)^2);
    yp_all(k+1) = term1 + theta_p(2)*up_all(k) + theta_p(3)*up_all(k-1);
    
    % 数字孪生线性方程
    % yd(t+1) = ...
    yd_all(k+1) = theta_d(1)*yd_all(k) + theta_d(2)*yd_all(k-1) + theta_d(3)*ud_all(k) + theta_d(4)*ud_all(k-1);
    
    % --- (6) 更新积分状态 (为下一时刻准备) ---
    % z(t+1) = z(t) + r(t) - y(t)
    % *关键*：这里按照您的要求，使用方波 r(t) 进行误差计算
    zp_all(k+1) = zp_all(k) + (r_all(k) - yp_all(k));
    zd_all(k+1) = zd_all(k) + (r_all(k) - yd_all(k));
    
end

%% 4. 数据截取与绘图
% 截取有效数据段：从 offset+1 开始，长度为 len_plot (即 201 个点)
% 这样就能保证与 time (0:200) 长度完全一致
idx_start = offset + 1;
idx_end   = offset + len_plot;

outputr  = r_all(idx_start:idx_end);   % 方波信号
outputyr = yr_all(idx_start:idx_end);  % 期望输出
outputyp = yp_all(idx_start:idx_end);  % 物理输出
outputyd = yd_all(idx_start:idx_end);  % 数字输出

% 检查维度 (调试用)
fprintf('Time length: %d\n', length(time));
fprintf('Data length: %d\n', length(outputyp));

% 绘图
figure('Color', 'w', 'Position', [100, 100, 700, 500]);
hold on; box on; grid on;

% 1. 绘制方波输入 r(t) - 灰色虚线
plot(time, outputr, 'color', [0.6 0.6 0.6], 'LineStyle', '--', 'LineWidth', 1, 'DisplayName', 'Square Wave Input r(t)');

% 2. 绘制期望输出 yr(t) - 黑色实线 (Figure 11 中的 Reference)
plot(time, outputyr, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Desired Output y_r(t)');

% 3. 绘制物理输出 yp(t) - 蓝色实线
plot(time, outputyp, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Real Plant y_p(t)');

% 4. 绘制数字输出 yd(t) - 红色点划线
plot(time, outputyd, 'r-.', 'LineWidth', 1.5, 'DisplayName', 'Digital Twin y_d(t)');

% 添加 t=100 标注
x_arrow = 100;
ylim([-0.2 1.4]);
line([100 100], [-0.2 1.4], 'Color', 'k', 'LineStyle', ':');
text(102, 1.3, 'Parameters Retuned', 'FontSize', 10);

% 图例设置
legend('Location', 'northeast');
xlabel('Time (steps)');
ylabel('Outputs of the system');
title('Case 3: Performance Retention');
xlim([0 200]);

set(gca, 'FontSize', 10, 'FontName', 'Times New Roman');

fprintf('绘图完成。\n');
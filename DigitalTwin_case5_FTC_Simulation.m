%% 数字孪生控制器容错控制复现 (包含 Fig.14 和 Fig.15)
% 对应论文中的 Case 5: Controller-fault tolerant control
% 复现无容错(Fig.14) 和 有容错(Fig.15) 两种情况

clc; clear; close all;

%% === 全局参数设置 ===
T_steps = 200;  % 仿真总步数
N = T_steps + 10; % 数组长度

% 1. 物理对象参数 (Real Plant)
theta_p = [0.5, 1.5, 0.3];

% 2. 数字孪生参数 (Digital Plant - Estimated)
theta_d = [1.3178, -0.3866, 1.4238, -1.2606];

% 3. 控制器参数 (PI gains)
K_dp = 0.75;  K_di = 0.035; % D-Controller
K_pp = 0.75;  K_pi = 0.035; % R-Controller
h_d = 1;      h_p = 1;      % 积分项增益

%% === 预计算参考输入 r(t) ===
r = zeros(1, N);
for k = 1:N
    t_val = k - 1; 
    % 方波信号逻辑
    if (t_val >= 0 && t_val < 50) || (t_val >= 100 && t_val < 150)
        r(k) = 1;
    else
        r(k) = 0;
    end
end

%% === 运行两次仿真 ===
% mode = 1: 无容错控制 (对应 Fig. 14)
% mode = 2: 有容错控制 (对应 Fig. 15)

for mode = 1:2
    
    % --- 初始化状态变量 ---
    yp = zeros(1, N); % 物理对象输出
    up = zeros(1, N); % 物理对象控制输入
    zp = zeros(1, N); % R-Controller 积分
    
    yd = zeros(1, N); % 数字孪生输出
    ud = zeros(1, N); % 数字孪生控制输入
    zd = zeros(1, N); % D-Controller 积分
    
    yr = zeros(1, N); % 期望参考模型输出 (Desired Output)
    
    % --- 仿真循环 ---
    for k = 3:T_steps
        
        % 1. 参考模型更新 (Reference Model dynamics)
        % y_r(t+2) = ... (基于当前 k 计算未来的 k+1)
        yr(k+1) = 1.7210*yr(k) - 0.7558*yr(k-1) + 0.0182*r(k) + 0.0166*r(k-1);
        
        % 2. 数字孪生系统 (Digital Twin Loop)
        % 始终正常运行，不受物理控制器故障影响
        
        % D-Controller 积分状态 z_d
        zd(k+1) = zd(k) + h_d * (r(k) - yd(k));
        
        % D-Controller 输出 u_d
        num_d = (K_dp - theta_d(1))*yd(k) - theta_d(2)*yd(k-1) + K_di*zd(k) - theta_d(4)*ud(k-1);
        ud(k) = num_d / theta_d(3);
        
        % 数字对象状态更新 y_d
        yd(k+1) = theta_d(1)*yd(k) + theta_d(2)*yd(k-1) + theta_d(3)*ud(k) + theta_d(4)*ud(k-1);
        
        % 3. 物理实体系统 (Physical Plant Loop)
        
        % R-Controller 积分状态 z_p
        zp(k+1) = zp(k) + h_p * (r(k) - yp(k));
        
        % R-Controller 正常计算出的输出 u_p_calc
        num_p = (K_pp - theta_d(1))*yp(k) - theta_d(2)*yp(k-1) + K_pi*zp(k) - theta_d(4)*up(k-1);
        u_p_calc = num_p / theta_d(3);
        
        % --- 故障注入与控制逻辑选择 ---
        time_step = k - 1;
        
        if time_step < 80
            % [正常阶段]：物理控制器正常工作
            up(k) = u_p_calc;
        else
            % [故障阶段]：t >= 80，R-Controller 失效
            if mode == 1
                % Case 1 (Fig. 14): 无容错
                % 模拟控制器失效输出随机噪声 (-0.5 到 0.5)
                up(k) = -0.5 + rand(); 
            else
                % Case 2 (Fig. 15): 有容错 (Controller-Fault Tolerant)
                % 使用数字控制器的输出直接驱动物理对象
                up(k) = ud(k); 
            end
        end
        
        % 物理对象状态更新 (非线性模型) y_p
        % yp(t+1) = (theta_p1*yp(t)*yp(t-1))/(1+yp(t-1)^2) + ...
        term_nl = (theta_p(1) * yp(k) * yp(k-1)) / (1 + yp(k-1)^2);
        yp(k+1) = term_nl + theta_p(2) * up(k) + theta_p(3) * up(k-1);
    end
    
    % --- 存储绘图数据 ---
    t_axis = 0:T_steps;
    data_r = r(1:T_steps+1);
    data_yr = yr(1:T_steps+1);
    data_yd = yd(1:T_steps+1);
    data_yp = yp(1:T_steps+1);
    
    % --- 绘图 ---
    figure(mode);
    set(gcf, 'Color', 'w', 'Position', [100 + (mode-1)*600, 200, 560, 420]); % 并排显示窗口
    hold on; box on;
    
    % 1. Reference (方波) - 黑色虚线
    h1 = plot(t_axis, data_r, 'k-', 'LineWidth', 1.0); 
    
    % 2. Output y_r (期望轨迹) - 黑色实线 (细)
    h2 = plot(t_axis, data_yr, 'k-', 'LineWidth', 1.0);
    
    % 3. Output y_d (数字孪生) - 蓝色实线
    h3 = plot(t_axis, data_yd, 'b-', 'LineWidth', 1.5);
    
    % 4. Output y_p (物理实体) - 红色实线
    h4 = plot(t_axis, data_yp, 'r-', 'LineWidth', 1.5);
    
    % 标注故障线
    xline(80, 'k-', {'R-controller fails'}, 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1);
    
    % 图表装饰
    xlabel('Time (steps)');
    ylabel('Outputs of the system');
    ylim([-1.2, 1.5]); % 调整Y轴范围以容纳振荡
    xlim([0, 200]);
    grid on;
    
    % 图例
    legend([h1, h2, h3, h4], ...
           {'Reference r(t)', 'Output y_r(t)', 'Output y_d(t)', 'Output y_p(t)'}, ...
           'Location', 'best');
       
    if mode == 1
        title('Fig. 14: Without Controller-Fault Tolerant Control');
        % 箭头标注震荡
        text(120, 0.5, 'System Unstable', 'Color', 'red', 'FontSize', 10);
    else
        title('Fig. 15: With Controller-Fault Tolerant Control');
        % 箭头标注容错成功
        text(100, 1.1, 'FTC Active', 'Color', 'blue', 'FontSize', 10);
    end
    
    hold off;
    
end
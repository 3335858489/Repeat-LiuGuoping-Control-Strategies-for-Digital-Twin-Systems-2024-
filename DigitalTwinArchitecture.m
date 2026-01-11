%% 数字孪生控制系统 - 第二步：自适应模型跟踪 (Adaptive Model Tracking)
% 对应文献内容：Section III, Eq (23)-(32)
% 修复说明：
% 1. 修复了 'zd' 字段缺失的报错。
% 2. 实现了 RLS (递归最小二乘) 算法，在线更新 Theta_d 和 Phi_d。
% 环境：MATLAB 2025a (兼容旧版本)

clear; clc; close all;

%% 1. 全局初始化
T_final = 500;           % 仿真 50 秒
dt = 0.1;               % 采样时间
t = 0:dt:T_final;
N = length(t);

% 参考输入 r(t)
r = ones(1, N); 
r(1:5) = 0;             % 前 0.5s 为 0

%% 2. 物理子系统初始化 (Physical Subsystem - Truth)
% 这是一个稳定的二阶系统，作为我们的"被控对象"
phys.xp = zeros(2, N);  % 真实状态
phys.yp = zeros(1, N);  % 真实输出
phys.up = zeros(1, N);  % 真实输入
phys.zp = zeros(1, N);  % 真实控制器积分状态

% PID 参数
phys.kp = [0.5; 0.1];   % [Kp; Ki]
phys.hp = [1.0; dt];    % 积分器参数

% 真实物理参数 (这就是我们要辨识的目标)
% A_real: 特征值 < 1，保证物理系统本身是稳定的
A_real = [1.0,  0.1; 
         -0.3, 0.9];
B_real = [0; 0.1];
C_real = [1.0, 0];

%% 3. 虚拟子系统初始化 (Virtual Subsystem - Adaptive Model)
virt.xd = zeros(2, N);
virt.yd = zeros(1, N);
virt.ud = zeros(1, N);
virt.zd = zeros(1, N);  % <--- 修复点：显式初始化 zd，避免报错
virt.kd = phys.kp;      % 初始控制器参数相同
virt.hd = phys.hp;

% --- 自适应算法参数初始化 (对应 Eq 28-32) ---
% 1. 状态方程参数 Theta_d (对应 Eq 24)
% 包含 A 和 B 的估计值。维度：[x(k); u(k)] -> x(k+1)
% 真实 Theta 应该是 [A_real, B_real]' (转置后便于计算)
virt.Theta_d = zeros(3, 2); % 3行(2个状态+1个输入), 2列(输出是2个状态)
virt.Theta_d(1,1) = 0.8;    % 故意给错误的初值，测试自适应能力
virt.Theta_d(2,2) = 0.8;

% 2. 输出方程参数 Phi_d (对应 Eq 25)
% 维度：x(k) -> y(k)
virt.Phi_d = zeros(2, 1);   % 2行(状态), 1列(输出)
virt.Phi_d(1) = 0.5;        % 故意给错误的初值

% 3. RLS 协方差矩阵 P (对应文献中的 T_theta 和 T_phi)
% 初始化为大对角阵，表示初始不确定性大
virt.P_theta = 1000 * eye(3); 
virt.P_phi   = 1000 * eye(2); 

% 4. 遗忘因子 lambda (对应 Eq 28)
lambda = 0.98; % 0.95~0.99 之间，越小收敛越快但对噪声越敏感

%% 4. 主仿真循环
fprintf('开始自适应仿真...\n');

for k = 2:N
    % 当前时间步的参考信号
    r_curr = r(k-1);
    
    % =============================================================
    % Part A: 物理子系统运行 (产生真实数据)
    % =============================================================
    
    % 1. R-Controller 计算
    [phys.up(k-1), phys.zp(k)] = Controller_Func(...
        phys.zp(k-1), phys.yp(k-1), r_curr, phys.kp, phys.hp);
    
    % 2. Real Plant 状态更新 (真实物理过程)
    phys.xp(:, k) = A_real * phys.xp(:, k-1) + B_real * phys.up(k-1);
    
    % 3. Real Plant 输出
    phys.yp(k) = C_real * phys.xp(:, k);
    
    % =============================================================
    % Part B: 自适应参数辨识 (RLS Algorithm) - 核心部分
    % =============================================================
    % 文献思想：利用 t-1 时刻的数据预测 t 时刻，与 t 时刻真实值比较，更新参数
    
    if k > 2
        % --- 1. 辨识状态方程参数 Theta_d (Eq 29, 30) ---
        % 构建回归向量 F(t-1) = [x_p(t-1); u_p(t-1)]
        F_vec = [phys.xp(:, k-1); phys.up(k-1)]; 
        
        % 目标值：当前真实的物理状态 x_p(t)
        Target_State = phys.xp(:, k)'; % 转置为行向量以便计算
        
        % 调用 RLS 函数更新 Theta_d
        [virt.Theta_d, virt.P_theta] = RLS_Update(...
            virt.Theta_d, virt.P_theta, F_vec, Target_State, lambda);
            
        % --- 2. 辨识输出方程参数 Phi_d (Eq 31, 32) ---
        % 构建回归向量 G(t) = x_p(t) (文献中用的是 x_d 或 x_p，为了跟踪真实值，通常用真实状态)
        G_vec = phys.xp(:, k);
        
        % 目标值：当前真实的物理输出 y_p(t)
        Target_Output = phys.yp(k);
        
        % 调用 RLS 函数更新 Phi_d
        [virt.Phi_d, virt.P_phi] = RLS_Update(...
            virt.Phi_d, virt.P_phi, G_vec, Target_Output, lambda);
    end
    
    % =============================================================
    % Part C: 虚拟子系统运行 (使用更新后的参数)
    % =============================================================
    
    % 1. 解析当前学习到的参数
    % Theta_d 是 [A, B]' -> 所以 A = Theta_d(1:2, :)', B = Theta_d(3, :)'
    A_hat = virt.Theta_d(1:2, :)';
    B_hat = virt.Theta_d(3, :)';
    C_hat = virt.Phi_d';
    
    % 2. D-Controller 计算 (孪生控制器)
    [virt.ud(k-1), virt.zd(k)] = Controller_Func(...
        virt.zd(k-1), virt.yd(k-1), r_curr, virt.kd, virt.hd);
        
    % 3. Digital Plant 状态更新 (使用辨识出的 A_hat, B_hat)
    % Eq (24): x_d(t) = Theta' * F... 实际上就是 x = Ax + Bu
    virt.xd(:, k) = A_hat * virt.xd(:, k-1) + B_hat * virt.ud(k-1);
    
    % 4. Digital Plant 输出 (使用辨识出的 C_hat)
    virt.yd(k) = C_hat * virt.xd(:, k);
    
end

%% 5. 结果可视化
figure('Color', 'w', 'Position', [100, 100, 800, 600]);

% 5.1 输出跟踪对比
subplot(3,1,1);
plot(t, phys.yp, 'b-', 'LineWidth', 2); hold on;
plot(t, virt.yd, 'r--', 'LineWidth', 2);
plot(t, r, 'k:', 'LineWidth', 1);
title('自适应模型跟踪效果 (Adaptive Model Tracking)');
legend('Real Plant y_p', 'Digital Twin y_d (Adaptive)', 'Reference');
ylabel('Output'); grid on;
% 放大局部细节
axes('Position',[.6 .75 .25 .15])
box on
indexOfInterest = (t < 15) & (t > 0);
plot(t(indexOfInterest), phys.yp(indexOfInterest), 'b-', 'LineWidth', 1.5); hold on;
plot(t(indexOfInterest), virt.yd(indexOfInterest), 'r--', 'LineWidth', 1.5);
title('前15秒自适应过程'); grid on;

% 5.2 参数收敛过程 (A矩阵的第一个元素)
subplot(3,1,2);
% 真实值 A_real(1,1) = 1.0
line([0 T_final], [1.0 1.0], 'Color', 'b', 'LineStyle', '--'); hold on;
% 提取 Theta_d 的历史数据需要稍微改一下代码结构，这里简化展示
% 我们展示误差
err_track = phys.yp - virt.yd;
plot(t, err_track, 'k');
title('跟踪误差 e(t) = y_p - y_d');
ylabel('Error'); grid on;

% 5.3 控制输入
subplot(3,1,3);
plot(t, phys.up, 'b'); hold on;
plot(t, virt.ud, 'r--');
title('控制输入对比');
legend('Real Control', 'Digital Control');
xlabel('Time (s)'); grid on;


%% 6. 辅助函数 (Local Functions)

% --- 控制器函数 ---
function [u, z_next] = Controller_Func(z_curr, y, r, k_params, h_params)
    e = r - y;
    Kp = k_params(1);
    Ki = k_params(2);
    
    % 积分更新: z(k+1) = z(k) + e * dt
    % h_params(2) 存储的是 dt
    z_next = z_curr + e * h_params(2);
    
    % 控制律: u = Kp*e + Ki*z
    u = Kp*e + Ki*z_next;
end

% --- RLS (递归最小二乘) 核心算法 ---
% 对应 Eq (29)-(32)
% Theta: 当前参数估计值 (列向量或矩阵)
% P:     当前协方差矩阵 (方阵)
% phi:   回归向量 (Regressor Vector)
% y:     观测值 (真实值)
% lambda: 遗忘因子
function [Theta_new, P_new] = RLS_Update(Theta_old, P_old, phi, y, lambda)
    % 1. 计算增益向量 K (Eq 29/31 中的这一项: P*phi / (lambda + phi'*P*phi))
    % 分母是一个标量
    denominator = lambda + phi' * P_old * phi;
    K = (P_old * phi) / denominator;
    
    % 2. 计算预测误差 (Prediction Error)
    % 注意：如果 y 是向量，Theta 是矩阵，这里的维度要匹配
    % y: [1 x m], Theta: [n x m], phi: [n x 1] -> phi'*Theta: [1 x m]
    error = y - phi' * Theta_old;
    
    % 3. 更新参数 Theta (Eq 29/31)
    Theta_new = Theta_old + K * error;
    
    % 4. 更新协方差矩阵 P (Eq 30/32)
    % P_new = (I - K * phi') * P_old / lambda
    % 这里的 I 是单位阵，维度与 P 相同
    I = eye(size(P_old));
    P_new = (I - K * phi') * P_old / lambda;
end
%% 清除工作区和命令窗口
clear; clc; % 清除工作区变量和命令窗口内容

%% 定义输入信号
t = linspace(-5, 5, 1000); % 在 -5 到 5 之间生成 1000 个时间点
x = sinc(5 * t); % 使用 sinc 函数（5*t）定义输入信号 x
x_trans = sinc(5 * t) + exp(-100i * t) .* sinc(5 * t) + exp(100i * t) .* sinc(5 * t); % 对 x 进行频移，生成新的信号 x_trans，包含原始频率和 ±100 的频移频率
N = length(x_trans); % 获取信号 x_trans 的长度

%% 对信号进行傅里叶变换
X = fftshift(fft(x_trans)); % 对 x_trans 进行傅里叶变换，并将零频分量移到频谱中心

%% 绘制输入信号和傅里叶变换结果
figure('Name', '输入信号和傅里叶变换', 'Position', [100, 100, 1000, 800]); % 创建新图形窗口，命名为“输入信号和傅里叶变换”

% 绘制输入信号 x 的时域波形
subplot(2, 2, 1);
plot(t, x);
title('原始输入信号');
xlabel('时间 (s)');
ylabel('幅值');
grid on;

% 绘制频移后的信号 x_trans 的时域波形
subplot(2, 2, 2);
plot(t, real(x_trans)); % 取实部进行绘制
title('频移后的输入信号');
xlabel('时间 (s)');
ylabel('幅值');
grid on;

% 定义频率向量
fs = 1 / (t(2) - t(1)); % 计算采样频率
f = linspace(-fs/2, fs/2, N); % 定义频率向量，用于频域显示

% 绘制原始信号 x 的频域幅度谱
subplot(2, 2, 3);
XX = fftshift(fft(x)) / N; % 对原始信号 x 进行傅里叶变换并归一化
plot(f, abs(XX));
title('原始信号的傅里叶变换');
xlabel('频率 (Hz)');
ylabel('幅值');
grid on;

% 绘制频移后信号 x_trans 的频域幅度谱
subplot(2, 2, 4);
plot(f, abs(X));
title('频移后信号的傅里叶变换');
xlabel('频率 (Hz)');
ylabel('幅值');
grid on;

%% 设计并应用低通滤波器
fc = 5; % 设置低通滤波器的截止频率为 5 Hz
order = 50; % 设置滤波器阶数为 50
b = fir1(order, fc / (fs / 2)); % 使用 fir1 函数生成低通滤波器的系数，将截止频率归一化到 Nyquist 频率（fs/2）
y = filter(b, 1, x_trans); % 对频移后的信号 x_trans 进行滤波，得到滤波后的信号 y
Y = fftshift(fft(y)); % 对滤波后的信号 y 进行傅里叶变换，将零频分量移到频谱中心

%% 绘制输入输出信号及其傅里叶变换结果
figure('Name', '输入输出信号及其傅里叶变换', 'Position', [100, 100, 1000, 800]);

% 绘制原始信号 x_trans 的时域波形
subplot(4, 1, 1);
plot(t, real(x_trans)); % 取实部进行绘制
title('原始信号');
xlabel('时间 (s)');
ylabel('幅值');
grid on;

% 绘制滤波后的信号 y 的时域波形
subplot(4, 1, 2);
plot(t, real(y)); % 取实部进行绘制
title('滤波后的信号');
xlabel('时间 (s)');
ylabel('幅值');
grid on;

% 绘制原始信号 x_trans 的频域幅度谱
subplot(4, 1, 3);
plot(f, abs(X));
title('原始信号的傅里叶变换');
xlabel('频率 (Hz)');
ylabel('幅值');
grid on;

% 绘制滤波后信号 y 的频域幅度谱
subplot(4, 1, 4);
plot(f, abs(Y));
title('滤波后信号的傅里叶变换');
xlabel('频率 (Hz)');
ylabel('幅值');
grid on;

%% 绘制滤波器的频率响应
figure('Name', '滤波器的频率响应', 'Position', [100, 100, 800, 600]);
freqz(b, 1, 512, fs); % 绘制滤波器的频率响应特性，指定采样频率
title('滤波器频率响应');
xlabel('频率 (Hz)');
ylabel('幅值 (dB)');
grid on;

%% 调制和解调过程
% 调制
fc = 20000; % 设置载波频率为 20000 Hz
c = sin(2 * pi * fc * t); % 生成正弦载波信号
m = real(y) .* c; % 对滤波后的信号 y 进行调制，与载波信号 c 相乘得到调制后的信号 m

% 绘制调制前后的信号及其频谱
figure('Name', '信号调制过程', 'Position', [100, 100, 1000, 800]);

% 绘制调制前的信号（滤波后的信号）的时域波形
subplot(2, 2, 1);
plot(t, real(y));
title('调制前的信号');
xlabel('时间 (s)');
ylabel('幅值');
grid on;

% 绘制调制前信号的频域幅度谱
subplot(2, 2, 2);
Y_mod = fftshift(fft(real(y)));
plot(f, abs(Y_mod));
title('调制前信号的频谱');
xlabel('频率 (Hz)');
ylabel('幅值');
grid on;

% 绘制调制后的信号 m 的时域波形
subplot(2, 2, 3);
plot(t, m);
title('调制后的信号');
xlabel('时间 (s)');
ylabel('幅值');
grid on;

% 绘制调制后信号的频域幅度谱
subplot(2, 2, 4);
M = fftshift(fft(m));
plot(f, abs(M));
title('调制后信号的频谱');
xlabel('频率 (Hz)');
ylabel('幅值');
grid on;

% 解调
y_demod = m .* sin(2 * pi * fc * t); % 对调制后的信号进行解调，与载波信号相乘
baseband_sig = y_demod .* sin(2 * pi * fc * t); % 进一步处理解调后的信号
envelope = abs(hilbert(baseband_sig)); % 计算解调后信号的包络

% 滤波
[b_filt, a_filt] = butter(6, 0.5, 'low'); % 设计一个六阶低通巴特沃斯滤波器，截止频率为 0.5（归一化频率）
z = filter(b_filt, a_filt, y_demod); % 对解调后的信号进行滤波，得到滤波后的信号 z

% 绘制解调解调前后的信号和频谱
figure('Name', '信号解调过程', 'Position', [100, 100, 1000, 800]);

% 绘调解调前的信号（调制后的信号 m）的时域波形
subplot(3, 2, 1);
plot(t, m);
title('解调前的信号');
xlabel('时间 (s)');
ylabel('幅值');
grid on;

% 绘调解调前信号的频域幅度谱
subplot(3, 2, 3);
plot(f, abs(M));
title('解调前信号的频谱');
xlabel('频率 (Hz)');
ylabel('幅值');
grid on;

% 绘调解调后的信号的时域波形
subplot(3, 2, 2);
plot(t, real(y_demod)); % 取实部进行绘制
title('解调后的信号');
xlabel('时间 (s)');
ylabel('幅值');
grid on;

% 绘调解调后信号的频域幅度谱
subplot(3, 2, 4);
Y_demod = fftshift(fft(y_demod));
plot(f, abs(Y_demod));
title('解调后信号的频谱');
xlabel('频率 (Hz)');
ylabel('幅值');
grid on;

% 绘制滤波后信号 z 的频域幅度谱
subplot(3, 2, 5);
Z = fftshift(fft(z));
plot(f, abs(Z));
title('滤波后信号的频谱');
xlabel('频率 (Hz)');
ylabel('幅值');
grid on;

%% 拉普拉斯变换
syms s t;
laplace_x = laplace(sinc(5*t), t, s); % 对 sinc(5*t) 进行拉普拉斯变换
figure('Name', '拉普拉斯变换');
ezplot(abs(laplace_x), [0, 100]); % 绘制拉普拉斯变换结果的幅度谱，限定复频率范围
title('拉普拉斯变换幅度谱');
xlabel('复频率');
ylabel('幅值');
grid on;

%% 信号采样与 Z 变换
fs_sample = 10; % 设置采样频率为 10 Hz
t_sample = -5 : 1 / fs_sample : 5; % 生成采样时间点
x_sampled = sinc(5 * t_sample); % 对原始信号进行采样，得到采样信号 x_sampled

% 绘制采样信号的时域波形
figure('Name', '采样信号', 'Position', [100, 100, 800, 600]);
stem(t_sample, x_sampled, 'filled'); % 使用填充的茎状图显示采样点
title('采样信号');
xlabel('时间 (s)');
ylabel('幅值');
grid on;

%% 计算 Z 变换并绘制幅度谱
% 定义离散时间信号 x[n]
n = -5 * fs_sample : 5 * fs_sample; % 假设采样点范围与原始时间范围一致
x_n = sinc(5 * (n / fs_sample)); % 采样后的信号

% 计算 Z 变换的频率响应
N_fft = 1000; % 设置 FFT 点数为 1000
X_z = fftshift(abs(fft(x_n, N_fft))); % 使用指定 FFT 点数计算幅度谱

% 定义频率向量
w = linspace(0, pi, N_fft); % 定义数字频率范围（归一化频率），长度与 FFT 点数一致

% 绘制 Z 变换结果的幅度谱
figure('Name', 'Z 变换', 'Position', [100, 100, 800, 600]);
plot(w / pi, X_z);
title('Z 变换幅度谱');
xlabel('归一化频率 (\times \pi)');
ylabel('幅值');
grid on;

%% 信号恢复
t_recovered = linspace(-5, 5, 1000); % 定义恢复的时间范围，与原始信号的时间范围相同，生成 1000 个时间点
x_recovered = interp1(t_sample, x_sampled, t_recovered, 'spline'); % 使用样条插值对采样信号进行恢复，得到恢复后的信号 x_recovered

% 绘制恢复后的信号的时域波形
figure('Name', '信号恢复', 'Position', [100, 100, 800, 600]);
plot(t_recovered, x_recovered); % 绘制恢复后的信号的时域波形
hold on;
stem(t_sample, x_sampled, 'r', 'filled'); % 在恢复后的信号上绘制原始采样点，使用红色填充茎状图显示
title('恢复后的信号');
xlabel('时间 (s)');
ylabel('幅值');
legend('恢复后的信号', '采样点', 'Location', 'best');
grid on;
hold off;

%% 希尔伯特变换分析
analytic_signal = hilbert(real(y)); % 对滤波后的信号进行希尔伯特变换
amplitude_envelope = abs(analytic_signal); % 计算信号的幅度包络
instantaneous_phase = unwrap(angle(analytic_signal)); % 计算信号的瞬时相位
t = linspace(-5, 5, 1000); % 定义时间范围，从 -5 到 5，生成 1000 个时间点
figure('Name', '希尔伯特变换分析', 'Position', [100, 100, 800, 600]);
subplot(3,1,1);
plot(t, real(y));
title('原始信号');
grid on;
subplot(3,1,2);
plot(t, amplitude_envelope);
title('信号幅度包络');
grid on;
subplot(3,1,3);
plot(t, instantaneous_phase);
title('信号瞬时相位');
grid on;

%% 综合分析
figure('Name', '综合分析', 'Position', [100, 100, 1000, 800]);

% 绘制原始信号、频移信号、滤波信号的对比
t = linspace(-5, 5, 1000); % 定义时间范围，从 -5 到 5，生成 1000 个时间点

subplot(3, 1, 1);
plot(t, real(x_trans)); % 取实部进行绘制
title('原始频移信号');
xlabel('时间 (s)');
ylabel('幅值');
grid on;

subplot(3, 1, 2);
plot(t, real(y)); % 取实部进行绘制
title('滤波后信号');
xlabel('时间 (s)');
ylabel('幅值');
grid on;

subplot(3, 1, 3);
plot(t_recovered, x_recovered);
title('恢复后的信号');
xlabel('时间 (s)');
ylabel('幅值');
grid on;

%% 保持图形窗口可见
drawnow;
%% 清除工作区和命令窗口
clear; clc;

%% 定义输入信号
t = linspace(-5, 5, 1000);
x = sinc(5 * t);
x_trans = sinc(5 * t) + exp(-100i * t) .* sinc(5 * t) + exp(100i * t) .* sinc(5 * t);
N = length(x_trans);

%% 对信号进行傅里叶变换
X = fftshift(fft(x_trans));

%% 绘制输入信号和傅里叶变换结果
fig1 = figure;
set(fig1, 'Name', '输入信号和傅里叶变换', 'Position', [100, 100, 1000, 800]);

subplot(2, 2, 1);
plot(t, x);
title('原始输入信号');
xlabel('时间 (s)');
ylabel('幅值');
grid on;

subplot(2, 2, 2);
plot(t, real(x_trans));
title('频移后的输入信号');
xlabel('时间 (s)');
ylabel('幅值');
grid on;

fs = 1 / (t(2) - t(1));
f = linspace(-fs/2, fs/2, N);

subplot(2, 2, 3);
XX = fftshift(fft(x)) / N;
plot(f, abs(XX));
title('原始信号的傅里叶变换');
xlabel('频率 (Hz)');
ylabel('幅值');
grid on;

subplot(2, 2, 4);
plot(f, abs(X));
title('频移后信号的傅里叶变换');
xlabel('频率 (Hz)');
ylabel('幅值');
grid on;

%% 设计并应用低通滤波器
fc = 5;
order = 50;
b = fir1(order, fc / (fs / 2));
y = filter(b, 1, x_trans);
Y = fftshift(fft(y));

%% 绘制输入输出信号及其傅里叶变换结果
fig2 = figure;
set(fig2, 'Name', '输入输出信号及其傅里叶变换', 'Position', [100, 100, 1000, 800]);

subplot(4, 1, 1);
plot(t, real(x_trans));
title('原始信号');
xlabel('时间 (s)');
ylabel('幅值');
grid on;

subplot(4, 1, 2);
plot(t, real(y));
title('滤波后的信号');
xlabel('时间 (s)');
ylabel('幅值');
grid on;

subplot(4, 1, 3);
plot(f, abs(X));
title('原始信号的傅里叶变换');
xlabel('频率 (Hz)');
ylabel('幅值');
grid on;

subplot(4, 1, 4);
plot(f, abs(Y));
title('滤波后信号的傅里叶变换');
xlabel('频率 (Hz)');
ylabel('幅值');
grid on;

%% 绘制滤波器的频率响应
fig3 = figure;
set(fig3, 'Name', '滤波器的频率响应', 'Position', [100, 100, 800, 600]);
freqz(b, 1, 512, fs);
title('滤波器频率响应');
xlabel('频率 (Hz)');
ylabel('幅值 (dB)');
grid on;

%% 调制和解调过程
fc = 20000;
c = sin(2 * pi * fc * t);
m = real(y) .* c;

fig4 = figure;
set(fig4, 'Name', '信号调制过程', 'Position', [100, 100, 1000, 800]);

subplot(2, 2, 1);
plot(t, real(y));
title('调制前的信号');
xlabel('时间 (s)');
ylabel('幅值');
grid on;

subplot(2, 2, 2);
Y_mod = fftshift(fft(real(y)));
plot(f, abs(Y_mod));
title('调制前信号的频谱');
xlabel('频率 (Hz)');
ylabel('幅值');
grid on;

subplot(2, 2, 3);
plot(t, m);
title('调制后的信号');
xlabel('时间 (s)');
ylabel('幅值');
grid on;

subplot(2, 2, 4);
M = fftshift(fft(m));
plot(f, abs(M));
title('调制后信号的频谱');
xlabel('频率 (Hz)');
ylabel('幅值');
grid on;

% 解调
y_demod = m .* sin(2 * pi * fc * t);
baseband_sig = y_demod .* sin(2 * pi * fc * t);
envelope = abs(hilbert(baseband_sig));

% 滤波
[b_filt, a_filt] = butter(6, 0.5, 'low');
z = filter(b_filt, a_filt, y_demod);

% 绘制解调解调前后的信号和频谱
fig5 = figure;
set(fig5, 'Name', '信号解调过程', 'Position', [100, 100, 1000, 800]);

subplot(3, 2, 1);
plot(t, m);
title('解调前的信号');
xlabel('时间 (s)');
ylabel('幅值');
grid on;

subplot(3, 2, 3);
plot(f, abs(M));
title('解调前信号的频谱');
xlabel('频率 (Hz)');
ylabel('幅值');
grid on;

subplot(3, 2, 2);
plot(t, real(y_demod));
title('解调后的信号');
xlabel('时间 (s)');
ylabel('幅值');
grid on;

subplot(3, 2, 4);
Y_demod = fftshift(fft(y_demod));
plot(f, abs(Y_demod));
title('解调后信号的频谱');
xlabel('频率 (Hz)');
ylabel('幅值');
grid on;

subplot(3, 2, 5);
Z = fftshift(fft(z));
plot(f, abs(Z));
title('滤波后信号的频谱');
xlabel('频率 (Hz)');
ylabel('幅值');
grid on;

%% 拉普拉斯变换 (修改部分)
syms s t_real;
% 显式定义符号化的 sinc 函数
sinc_sym = sin(5*pi*t_real) / (5*pi*t_real); 
laplace_x = laplace(sinc_sym, t_real, s);

fig6 = figure;
set(fig6, 'Name', '拉普拉斯变换', 'Position', [100, 100, 800, 600]);
ezplot(abs(laplace_x), [0, 100]); 
title('拉普拉斯变换幅度谱');
xlabel('复频率');
ylabel('幅值');
grid on;

%% 信号采样与 Z 变换
fs_sample = 10;
t_sample = -5 : 1 / fs_sample : 5;
x_sampled = sinc(5 * t_sample);

% 绘制采样信号的时域波形
fig7 = figure;
set(fig7, 'Name', '采样信号', 'Position', [100, 100, 800, 600]);
stem(t_sample, x_sampled, 'MarkerFaceColor', 'auto');
title('采样信号');
xlabel('时间 (s)');
ylabel('幅值');
grid on;

%% 计算 Z 变换并绘制幅度谱
n = -5 * fs_sample : 5 * fs_sample;
x_n = sinc(5 * (n / fs_sample));
N_fft = 1000;
X_z = fftshift(abs(fft(x_n, N_fft)));
w = linspace(0, pi, N_fft);

fig8 = figure;
set(fig8, 'Name', 'Z 变换', 'Position', [100, 100, 800, 600]);
plot(w / pi, X_z);
title('Z 变换幅度谱');
xlabel('归一化频率 (\times \pi)');
ylabel('幅值');
grid on;

%% 信号恢复
t_recovered = linspace(-5, 5, 1000);
x_recovered = interp1(t_sample, x_sampled, t_recovered, 'spline');

fig9 = figure;
set(fig9, 'Name', '信号恢复', 'Position', [100, 100, 800, 600]);
plot(t_recovered, x_recovered);
hold on;
stem(t_sample, x_sampled, 'MarkerFaceColor', 'r');
title('恢复后的信号');
xlabel('时间 (s)');
ylabel('幅值');
legend('恢复后的信号', '采样点', 'Location', 'best');
grid on;
hold off;

%% 希尔伯特变换分析 (移除冲突的 t 定义)
analytic_signal = hilbert(real(y));
amplitude_envelope = abs(analytic_signal);
instantaneous_phase = unwrap(angle(analytic_signal));

fig10 = figure;
set(fig10, 'Name', '希尔伯特变换分析', 'Position', [100, 100, 800, 600]);
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
fig11 = figure;
set(fig11, 'Name', '综合分析', 'Position', [100, 100, 1000, 800]);

subplot(3, 1, 1);
plot(t, real(x_trans));
title('原始频移信号');
xlabel('时间 (s)');
ylabel('幅值');
grid on;

subplot(3, 1, 2);
plot(t, real(y));
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
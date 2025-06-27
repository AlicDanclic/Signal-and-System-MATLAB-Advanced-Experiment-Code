clear; clc;

%% ================ 所有变量计算 ================
% 基础信号
t = linspace(-5, 5, 1000);
x = sinc(5 * t);
xt = x + exp(-100i * t).*x + exp(100i * t).*x;
N = length(xt);
fs = 1 / (t(2) - t(1));
f = linspace(-fs/2, fs/2, N);

% 傅里叶变换
Xf = fftshift(fft(xt));
XX = fftshift(fft(x)) / N;

% 滤波器处理
fc = 5; ord = 50;
b = fir1(ord, fc / (fs / 2));
y = filter(b, 1, xt);
Yf = fftshift(fft(y));

% 调制处理
fc_mod = 20000;
c = sin(2 * pi * fc_mod * t);
m = real(y) .* c;
Mf = fftshift(fft(m));
Ymod = fftshift(fft(real(y)));

% 解调处理
yd = m .* sin(2 * pi * fc_mod * t);
[b_filt, a_filt] = butter(6, 0.5, 'low');
z = filter(b_filt, a_filt, yd);
Ydemod = fftshift(fft(yd));
Zf = fftshift(fft(z));

% LT变换
syms s t_real;
sinc_sym = sin(5*pi*t_real)/(5*pi*t_real); 
Lx = laplace(sinc_sym, t_real, s);

% 采样处理
fsamp = 10; 
tsamp = -5 : 1/fsamp : 5;
xsamp = sinc(5 * tsamp);

% Z变换分析
n = -5*fsamp : 5*fsamp;
xn = sinc(5*(n/fsamp));
Nfft = 1000;
Xz = fftshift(abs(fft(xn, Nfft)));
w = linspace(0, pi, Nfft);

% 信号重建
trec = linspace(-5, 5, 1000);
xrec = interp1(tsamp, xsamp, trec, 'spline');

% Hilbert分析
asig = hilbert(real(y));
env = abs(asig);
iphase = unwrap(angle(asig));

%% ================ 绘图部分 ================
%% 图1：信号分析
figure('Name','信号分析','Position',[100,100,1000,800]);

subplot(2,2,1);
plot(t, x); title('输入信号'); xlabel('t(s)'); ylabel('幅值'); grid on;

subplot(2,2,2);
plot(t, real(xt)); title('频移信号'); xlabel('t(s)'); ylabel('幅值'); grid on;

subplot(2,2,3);
plot(f, abs(XX)); title('输入频谱'); xlabel('f(Hz)'); ylabel('幅值'); grid on;

subplot(2,2,4);
plot(f, abs(Xf)); title('频移频谱'); xlabel('f(Hz)'); ylabel('幅值'); grid on;

%% 图2：滤波分析
figure('Name','滤波分析','Position',[100,100,1000,800]);

subplot(4,1,1);
plot(t, real(xt)); title('原始信号'); xlabel('t(s)'); ylabel('幅值'); grid on;

subplot(4,1,2);
plot(t, real(y)); title('滤波信号'); xlabel('t(s)'); ylabel('幅值'); grid on;

subplot(4,1,3);
plot(f, abs(Xf)); title('原始频谱'); xlabel('f(Hz)'); ylabel('幅值'); grid on;

subplot(4,1,4);
plot(f, abs(Yf)); title('滤波频谱'); xlabel('f(Hz)'); ylabel('幅值'); grid on;

%% 图3：滤波器响应
figure('Name','滤波器响应','Position',[100,100,800,600]);
freqz(b, 1, 512, fs); title('滤波器响应'); xlabel('f(Hz)'); ylabel('幅值(dB)'); grid on;

%% 图4：调制分析
figure('Name','调制分析','Position',[100,100,1000,800]);

subplot(2,2,1);
plot(t, real(y)); title('基带信号'); xlabel('t(s)'); ylabel('幅值'); grid on;

subplot(2,2,2);
plot(f, abs(Ymod)); title('基带频谱'); xlabel('f(Hz)'); ylabel('幅值'); grid on;

subplot(2,2,3);
plot(t, m); title('调制信号'); xlabel('t(s)'); ylabel('幅值'); grid on;

subplot(2,2,4);
plot(f, abs(Mf)); title('调制频谱'); xlabel('f(Hz)'); ylabel('幅值'); grid on;

%% 图5：解调分析
figure('Name','解调分析','Position',[100,100,1000,800]);

subplot(3,2,1);
plot(t, m); title('调制信号'); xlabel('t(s)'); ylabel('幅值'); grid on;

subplot(3,2,3);
plot(f, abs(Mf)); title('调制频谱'); xlabel('f(Hz)'); ylabel('幅值'); grid on;

subplot(3,2,2);
plot(t, real(yd)); title('解调信号'); xlabel('t(s)'); ylabel('幅值'); grid on;

subplot(3,2,4);
plot(f, abs(Ydemod)); title('解调频谱'); xlabel('f(Hz)'); ylabel('幅值'); grid on;

subplot(3,2,5);
plot(f, abs(Zf)); title('滤波频谱'); xlabel('f(Hz)'); ylabel('幅值'); grid on;

%% 图6：LT变换
figure('Name','LT变换','Position',[100,100,800,600]);
ezplot(abs(Lx), [0,100]); title('LT谱'); xlabel('s域'); ylabel('幅值'); grid on;

%% 图7：采样信号
figure('Name','采样信号','Position',[100,100,800,600]);
stem(tsamp, xsamp, 'MarkerFaceColor','auto'); title('采样信号'); xlabel('t(s)'); ylabel('幅值'); grid on;

%% 图8：Z变换
figure('Name','Z变换','Position',[100,100,800,600]);
plot(w/pi, Xz); title('Z谱'); xlabel('归一化频率'); ylabel('幅值'); grid on;

%% 图9：信号重建
figure('Name','信号重建','Position',[100,100,800,600]);
plot(trec, xrec); hold on;
stem(tsamp, xsamp, 'MarkerFaceColor','r');
title('重建信号'); xlabel('t(s)'); ylabel('幅值');
legend('重建','采样点'); grid on; hold off;

%% 图10：Hilbert分析
figure('Name','Hilbert分析','Position',[100,100,800,600]);

subplot(3,1,1);
plot(t, real(y)); title('原始信号'); grid on;

subplot(3,1,2);
plot(t, env); title('包络'); grid on;

subplot(3,1,3);
plot(t, iphase); title('瞬时相位'); grid on;

%% 图11：综合分析
figure('Name','综合分析','Position',[100,100,1000,800]);

subplot(3,1,1);
plot(t, real(xt)); title('原始频移信号'); xlabel('t(s)'); ylabel('幅值'); grid on;

subplot(3,1,2);
plot(t, real(y)); title('滤波后信号'); xlabel('t(s)'); ylabel('幅值'); grid on;

subplot(3,1,3);
plot(trec, xrec); title('恢复后的信号'); xlabel('t(s)'); ylabel('幅值'); grid on;

drawnow;
clc; clear all; close all;

%% Parameters
fs      = 500;          % sampling frequency (Hz)
T       = 10;           % duration (s)
N       = fs * T;       % number of samples
f       = 1.2;          % average heart rate (Hz) ~ 72 bpm

t = (0:N-1)/fs;         % time axis (s)

%% Define Gaussian model parameters for TWO beat types: normal (state 1) and abnormal (state 2)
% Each of alpha, b, theta is a struct array with fields x,y,z

% -------- State 1: Normal beat template --------
alpha(1).x = [0.2,  1.0,  0.3];    % P, QRS, T amplitudes (x)
alpha(1).y = [0.1, -0.5,  0.2];    % (y)
alpha(1).z = [0.05, 0.8,  0.25];   % (z)

b(1).x = [0.15, 0.05, 0.12];       % widths in phase (x)
b(1).y = [0.15, 0.05, 0.12];       % (y)
b(1).z = [0.15, 0.05, 0.12];       % (z)

theta(1).x = [-2.0, -0.1, 1.5];    % phase centers (x)
theta(1).y = [-2.0, -0.1, 1.5];    % (y)
theta(1).z = [-2.0, -0.1, 1.5];    % (z)

% -------- State 2: Abnormal beat template (e.g., PVC-like) --------
alpha(2).x = [0.1,  1.5,  0.1];    % larger, sharper QRS, smaller P/T
alpha(2).y = [0.05, -0.8, 0.1];
alpha(2).z = [0.02,  1.2, 0.05];

b(2).x = [0.12, 0.03, 0.10];       % narrower QRS
b(2).y = [0.12, 0.03, 0.10];
b(2).z = [0.12, 0.03, 0.10];

theta(2).x = [-2.0, 0.0, 1.6];     % slightly shifted timing
theta(2).y = [-2.0, 0.0, 1.6];
theta(2).z = [-2.0, 0.0, 1.6];

theta0 = -pi;                      % initial phase

%% Markov state transition matrix (STM) for beat types
% Rows sum to 1, S_ij = P(next state = j | current state = i)
% State 1: normal, State 2: abnormal
STM = [0.95 0.05;    % from normal: stay normal 95%, go abnormal 5%
       0.60 0.40];   % from abnormal: return to normal 60%, stay abnormal 40%

S0 = 2;              % initial state: normal beat

%% Generate abnormal 3D dipole / VCG with Markov-switching morphology
[vcg, phi] = vcg_gen_abnormal(N, fs, f, alpha, b, theta, theta0, STM, S0);

%% Define a lead-field matrix to get multichannel ECG from VCG
% L maps [x;y;z] -> ECG leads
L = [  0.6   0.1   0.2;   % Lead 1
      -0.3   0.8   0.1;   % Lead 2
       0.1   0.5  -0.4;   % Lead 3
       0.4  -0.2   0.7;   % Lead 4
      -0.5  -0.4   0.2;   % Lead 5
       0.2   0.3   0.6;   % Lead 6
       0.1  -0.7   0.3;   % Lead 7
      -0.2   0.2   0.5];  % Lead 8

% Form 3 x N matrix of dipole
D = [vcg.x; vcg.y; vcg.z];     % size: 3 x N

% Multichannel ECG: (nLeads x N)
ecg = L * D;

%% ---------- Version 1: Add small Gaussian measurement noise ----------
noise_level   = 0.02;          % adjust as needed
ecg_gaussian  = ecg + noise_level * randn(size(ecg));

[nLeads, ~] = size(ecg_gaussian);

figure; hold on;
offset = 1.5 * max(abs(ecg_gaussian(:)));  % vertical spacing

for k = 1:nLeads
    plot(t, ecg_gaussian(k,:) + (k-1)*offset);
end

xlabel('Time (s)');
ylabel('Amplitude + offset');
title('Abnormal multichannel ECG (Markov VCG) with Gaussian noise');
set(gca, 'YTick', (0:nLeads-1)*offset, ...
         'YTickLabel', arrayfun(@(k)sprintf('Lead %d',k),1:nLeads,'UniformOutput',false));
grid on;
saveas(gcf, 'abnormal_gaussian_noise.png');

%% ---------- Version 2: Add realistic muscle artifact (biosignal\_noise\_gen, mode 2) ----------
snr_ma = 2;                     % SNR in dB (signal / muscle artifact)
ecg_ma = zeros(size(ecg));

for k = 1:nLeads
    SignalPower = mean(ecg(k,:).^2);
    ma_noise    = biosignal_noise_gen(2, SignalPower, snr_ma, N, fs);  % real MA
    ecg_ma(k,:) = ecg(k,:) + ma_noise.';
end

figure; hold on;
offset = 1.5 * max(abs(ecg_ma(:)));

for k = 1:nLeads
    plot(t, ecg_ma(k,:) + (k-1)*offset);
end

xlabel('Time (s)');
ylabel('Amplitude + offset');
title(sprintf('Abnormal multichannel ECG with real muscle artifact (SNR = %g dB)', snr_ma));
set(gca, 'YTick', (0:nLeads-1)*offset, ...
         'YTickLabel', arrayfun(@(k)sprintf('Lead %d',k),1:nLeads,'UniformOutput',false));
grid on;
saveas(gcf, 'abnormal_muscle_noise.png');

%% ---------- Version 3: Add baseline wander (biosignal\_noise\_gen, mode 4) ----------
snr_bw = 0;                     % SNR in dB (ECG / baseline wander)

SignalPower_global = mean(ecg(:).^2);
bw_noise = biosignal_noise_gen(4, SignalPower_global, snr_bw, N, fs);   % single BW trace

ecg_bw = zeros(size(ecg));
for k = 1:nLeads
    ecg_bw(k,:) = ecg(k,:) + bw_noise.';      % same baseline drift on all leads
end

figure; hold on;
offset = 1.5 * max(abs(ecg_bw(:)));

for k = 1:nLeads
    plot(t, ecg_bw(k,:) + (k-1)*offset);
end

xlabel('Time (s)');
ylabel('Amplitude + offset');
title(sprintf('Abnormal multichannel ECG with baseline wander (SNR = %g dB)', snr_bw));
set(gca, 'YTick', (0:nLeads-1)*offset, ...
         'YTickLabel', arrayfun(@(k)sprintf('Lead %d',k),1:nLeads,'UniformOutput',false));
grid on;
saveas(gcf, 'abnormal_baselinewander_noise.png');

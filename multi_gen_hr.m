clc; clear all; close all;

%% Parameters
fs       = 500;          % sampling frequency (Hz)
baseHR   = 72;           % nominal heart rate (bpm)
HRamp    = 10;           % HR modulation amplitude (bpm)
nBeats   = 15;           % number of beats to simulate

% Heart rate vector per beat (slow sinusoidal variation around baseHR)
beatIdx  = 0:nBeats-1;
HR       = baseHR + HRamp * sin(2*pi*beatIdx/nBeats);   % 1 x nBeats (bpm)

%% Gaussian model parameters (simple P-QRS-T template)

% Amplitudes for x,y,z (three Gaussians: P, QRS, T)
alpha.x = [0.2,  1.0,  0.3];
alpha.y = [0.1, -0.5,  0.2];
alpha.z = [0.05, 0.8,  0.25];

% Widths of Gaussians (in radians of phase)
b.x = [0.15, 0.05, 0.12];
b.y = [0.15, 0.05, 0.12];
b.z = [0.15, 0.05, 0.12];

% Centers of Gaussians in phase (P, QRS, T positions)
theta.x = [-2.0, -0.1, 1.5];
theta.y = [-2.0, -0.1, 1.5];
theta.z = [-2.0, -0.1, 1.5];

theta0 = -pi;            % initial phase

%% Generate 3D dipole / VCG with variable heart rate
% Use MATRIX implementation (default, but we pass it explicitly)
[vcg, phi] = vcg_gen_var_hr(HR, fs, alpha, b, theta, theta0, 'MATRIX');

% Time axis based on length of generated signal
N = length(vcg.x);
t = (0:N-1)/fs;

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
D = [vcg.x; vcg.y; vcg.z];   % size: 3 x N

% Multichannel ECG: (nLeads x N)
ecg = L * D;
[nLeads, ~] = size(ecg);

%% ---------- Version 1: Gaussian measurement noise ----------
noise_level   = 0.02;               % adjust as needed
ecg_gaussian  = ecg + noise_level * randn(size(ecg));

figure; hold on;
offset = 1.5 * max(abs(ecg_gaussian(:)));  % vertical spacing

for k = 1:nLeads
    plot(t, ecg_gaussian(k,:) + (k-1)*offset);
end

xlabel('Time (s)');
ylabel('Amplitude + offset');
title('Variable-HR multichannel ECG (VCG model) with Gaussian noise');
set(gca, 'YTick', (0:nLeads-1)*offset, ...
         'YTickLabel', arrayfun(@(k)sprintf('Lead %d',k),1:nLeads,'UniformOutput',false));
grid on;
saveas(gcf, 'varHR_gaussian_noise.png');

%% ---------- Version 2: Real muscle artifact (biosignal_noise_gen, mode 2) ----------
snr_ma = 2;                        % SNR in dB (signal / muscle artifact)
ecg_ma = zeros(size(ecg));

for k = 1:nLeads
    SignalPower = mean(ecg(k,:).^2);
    ma_noise    = biosignal_noise_gen(2, SignalPower, snr_ma, N, fs);  % mode 2: MA
    ecg_ma(k,:) = ecg(k,:) + ma_noise.';
end

figure; hold on;
offset = 1.5 * max(abs(ecg_ma(:)));

for k = 1:nLeads
    plot(t, ecg_ma(k,:) + (k-1)*offset);
end

xlabel('Time (s)');
ylabel('Amplitude + offset');
title(sprintf('Variable-HR multichannel ECG with real muscle artifact (SNR = %g dB)', snr_ma));
set(gca, 'YTick', (0:nLeads-1)*offset, ...
         'YTickLabel', arrayfun(@(k)sprintf('Lead %d',k),1:nLeads,'UniformOutput',false));
grid on;
saveas(gcf, 'varHR_muscle_noise.png');

%% ---------- Version 3: Baseline wander (biosignal_noise_gen, mode 4) ----------
snr_bw = 0;                        % SNR in dB (ECG / baseline wander)

SignalPower_global = mean(ecg(:).^2);
bw_noise = biosignal_noise_gen(4, SignalPower_global, snr_bw, N, fs);   % single BW trace

ecg_bw = zeros(size(ecg));
for k = 1:nLeads
    ecg_bw(k,:) = ecg(k,:) + bw_noise.';        % same drift on all leads
end

figure; hold on;
offset = 1.5 * max(abs(ecg_bw(:)));

for k = 1:nLeads
    plot(t, ecg_bw(k,:) + (k-1)*offset);
end

xlabel('Time (s)');
ylabel('Amplitude + offset');
title(sprintf('Variable-HR multichannel ECG with baseline wander (SNR = %g dB)', snr_bw));
set(gca, 'YTick', (0:nLeads-1)*offset, ...
         'YTickLabel', arrayfun(@(k)sprintf('Lead %d',k),1:nLeads,'UniformOutput',false));
grid on;
saveas(gcf, 'varHR_baselinewander_noise.png');

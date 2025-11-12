clc; clear all ; close all 
%% Parameters
fs      = 500;          % sampling frequency (Hz)
T       = 10;           % duration (s)
N       = fs * T;       % number of samples
f       = 1.2;          % average heart rate (Hz) ~ 72 bpm
f_dev   = 0.1;          % 10% HR variability

% --- Gaussian model parameters (very simple P-QRS-T example) ---

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

% Beat-to-beat stochastic deviations (percent)
delta_alpha  = 0.15;    % 15% amplitude variation
delta_b      = 0.10;    % 10% width variation
delta_theta  = 0.05;    % 5% phase-center variation

teta0 = -pi;            % initial phase

%% Generate stochastic 3D dipole / VCG
[vcg, phi] = vcg_gen_stochastic(N, fs, f, f_dev, ...
                                alpha, delta_alpha, ...
                                b,     delta_b, ...
                                theta, delta_theta, ...
                                teta0);

t = (0:N-1)/fs;         % time axis (s)

%% Define a lead-field matrix to get multichannel ECG from VCG
% L maps [x;y;z] -> ECG leads
% Here we just pick a simple 8-lead example. In practice, use a realistic lead-field.
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

%% Optional: add small measurement noise
noise_level = 0.02;          % adjust as needed
ecg_guassian = ecg + noise_level * randn(size(ecg));

%% Plot multichannel ECG (stacked)
[nLeads, ~] = size(ecg_guassian);
figure; hold on;
offset = 1.5 * max(abs(ecg_guassian(:)));  % vertical spacing

for k = 1:nLeads
    plot(t, ecg_guassian(k,:) + (k-1)*offset);
end

xlabel('Time (s)');
ylabel('Amplitude + offset');
title('Synthetic stochastic multichannel ECG from VCG model');
set(gca, 'YTick', (0:nLeads-1)*offset, ...
         'YTickLabel', arrayfun(@(k)sprintf('Lead %d',k),1:nLeads,'UniformOutput',false));
grid on;
saveas(gcf, 'multi_stochastistic_guassian_noise.png');

%%
%% Add realistic muscle artifact using biosignal_noise_gen

snr_ma = 2;   % Desired SNR in dB (signal / muscle artifact) – adjust as needed
ecg_ma = zeros(size(ecg));   % preallocate

for k = 1:nLeads
    % Signal power of this lead
    SignalPower = mean(ecg(k,:).^2);
    
    % Generate real muscle artifact noise (mode 2)
    % Usage: MA = biosignal_noise_gen(2, SignalPower, snr, N, fs, seed_optional);
    ma_noise = biosignal_noise_gen(2, SignalPower, snr_ma, N, fs);
    
    % biosignal_noise_gen returns a column vector; transpose to row
    ecg_ma(k,:) = ecg(k,:) + ma_noise.';
end


%%Plot multichannel ECG with muscle artifact
figure; hold on;
offset = 1.5 * max(abs(ecg_ma(:)));  % vertical spacing

for k = 1:nLeads
    plot(t, ecg_ma(k,:) + (k-1)*offset);
end

xlabel('Time (s)');
ylabel('Amplitude + offset');
title(sprintf('Synthetic multichannel ECG with real muscle artifact (SNR = %g dB)', snr_ma));
set(gca, 'YTick', (0:nLeads-1)*offset, ...
         'YTickLabel', arrayfun(@(k)sprintf('Lead %d',k),1:nLeads,'UniformOutput',false));
grid on;
saveas(gcf, 'multi_stochastistic_muscle_noise.png');

%%
%% Add baseline wander noise using biosignal_noise_gen (mode 4)

snr_bw = 0;  % Desired SNR in dB (ECG / baseline wander) – adjust as needed

% Use overall signal power across all leads to scale the baseline
SignalPower_global = mean(ecg(:).^2);

% Generate a single baseline wander noise trace
% Usage: BW = biosignal_noise_gen(4, SignalPower, snr, N, fs, seed_optional);
bw_noise = biosignal_noise_gen(4, SignalPower_global, snr_bw, N, fs);   % column vector

% Build baseline-wander-corrupted ECG
ecg_bw = zeros(size(ecg));
for k = 1:nLeads
    ecg_bw(k,:) = ecg(k,:) + bw_noise.';  % same baseline wander added to all leads
end

% Plot multichannel ECG with baseline wander

figure; hold on;
offset = 1.5 * max(abs(ecg_bw(:)));  % vertical spacing

for k = 1:nLeads
    plot(t, ecg_bw(k,:) + (k-1)*offset);
end

xlabel('Time (s)');
ylabel('Amplitude + offset');
title(sprintf('Synthetic multichannel ECG with baseline wander (SNR = %g dB)', snr_bw));
set(gca, 'YTick', (0:nLeads-1)*offset, ...
         'YTickLabel', arrayfun(@(k)sprintf('Lead %d',k),1:nLeads,'UniformOutput',false));
grid on;
saveas(gcf, 'multi_stochastistic_baselinewander_noise.png');



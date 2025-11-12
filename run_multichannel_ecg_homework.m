clc ;clear ;close all 
run_multichannel_ecg_homework()

function run_multichannel_ecg_homework()
% Multichannel ECG Modeling and Stochastic Extensions (OSET-based)
% - 3D VCG maternal + fetal with stochastic variability
% - Projection to 12-lead ECG
% - Add realistic noise
% - Compute metrics (SNR, morphology variability, HRV, cross-lead corr)
% - Export figures and a Markdown report
%
% Requirements: Your OSET functions on the MATLAB path:
%   vcg_gen_stochastic, vcg_gen_state_space (optional), etc.

%% 0) Reproducibility + output dir
rng(7);
outdir = "outputs";
if ~exist(outdir, "dir"); mkdir(outdir); end

%% 1) Global parameters
fs = 500;                % Hz
dur_s = 30;              % seconds
N = dur_s*fs;

% --- Maternal parameters (stochastic VCG) ---
f_m_bpm      = 70;       % mean HR (bpm)
fdev_m       = 0.06;     % +/- 6% beat-to-beat HR jitter
theta0_m     = -2.2;     % initial phase
d_alpha_m    = 0.10;     % per-beat jitter (%)
d_b_m        = 0.10;
d_theta_m    = 0.05;
alpha_m = struct('x', [0.8 0.5 0.3], 'y', [0.6 0.4 0.2], 'z', [0.7 0.5 0.25]);
b_m     = struct('x', [0.06 0.05 0.05], 'y', [0.06 0.05 0.05], 'z', [0.06 0.05 0.05]);
theta_m = struct('x', [-1.10 -0.15 0.20], 'y', [-1.05 -0.20 0.22], 'z', [-1.00 -0.25 0.24]);

% --- Fetal parameters (stochastic VCG) ---
f_f_bpm   = 150;
fdev_f    = 0.08;
theta0_f  = 0.9;
d_alpha_f = 0.12;
d_b_f     = 0.12;
d_theta_f = 0.06;
alpha_f = struct('x', 0.35*[1 0.7 0.5], 'y', 0.30*[1 0.7 0.5], 'z', 0.25*[1 0.7 0.5]);
b_f     = struct('x', [0.05 0.045 0.04], 'y', [0.05 0.045 0.04], 'z', [0.05 0.045 0.04]);
theta_f = struct('x', [-1.00 -0.10 0.18], 'y', [-0.95 -0.16 0.20], 'z', [-0.90 -0.20 0.22]);

% --- Multichannel projection (lead field) ---
% For realism, substitute an inverse-Dower 12x3 if you have it.
A = randn(12,3); [Q,~] = qr(A,0); L_m = Q;            % maternal lead field (12x3)
R = axang2rotm([0 0 1 deg2rad(20)]);                  % slight rotation for fetal orientation
L_f = Q*R*diag([0.8 0.8 0.8]);                        % scaled/rotated fetal field

% --- Noise controls ---
target_snr_db = 15;          % overall SNR target
bw_amp = 0.05;               % baseline wander scale (~5%)
emg_scale = 0.02;            % EMG-like noise scale

%% 2) Synthesize maternal & fetal VCGs with stochastic variability
[vcg_m, phi_m] = vcg_gen_stochastic( ...
    N, fs, f_m_bpm/60, fdev_m, alpha_m, d_alpha_m, b_m, d_b_m, theta_m, d_theta_m, theta0_m);

[vcg_f, phi_f] = vcg_gen_stochastic( ...
    N, fs, f_f_bpm/60, fdev_f, alpha_f, d_alpha_f, b_f, d_b_f, theta_f, d_theta_f, theta0_f);

%% 3) Project to 12-lead and mix maternal + fetal
ecg_m = project_vcg_to_ecg(vcg_m, L_m);            % 12 x N
ecg_f = project_vcg_to_ecg(vcg_f, L_f);
mix_clean = ecg_m + 0.35*ecg_f;                    % fetal ~35% amplitude of maternal

%% 4) Add realistic composite noise
ecg_noisy = add_ecg_noise(mix_clean, fs, target_snr_db, bw_amp, emg_scale);

%% 5) Compute metrics
% SNR per lead
[snr_per_lead_db, mean_snr_db] = snr_per_lead(mix_clean, ecg_noisy);

% Choose a maternal-dominant lead (max variance proxy)
[~, lead_m_dom] = max(var(ecg_m,0,2));
x_m = ecg_noisy(lead_m_dom,:);

% Simple R-peak detection for maternal (bandpass + threshold)
rpos_m = simple_rpeaks(x_m, fs);

% Morphology variability on maternal-dominant lead (cosine distance between consecutive beats)
morph_var = beat_to_beat_variability(x_m, rpos_m, fs);

% Maternal HRV (SDNN)
rr_m_ms = diff(rpos_m)/fs*1000;
sdnn_m = std(rr_m_ms);

% Very crude fetal R-peak proxy: high-pass + detect around 100–250 bpm
x_f_hp = bandpass(ecg_noisy(lead_m_dom,:), [40 90], fs);  % emphasize fetal content
rpos_f = simple_rpeaks(x_f_hp, fs, 'expected_bpm', [110 200]);
rr_f_ms = diff(rpos_f)/fs*1000;
sdnn_f = std(rr_f_ms);

% Cross-lead correlation (average of upper triangle)
clcMat = corrcoef(ecg_noisy.'); % leads as variables
clc_upper = clcMat(triu(true(size(clcMat)),1));
crosslead_corr = mean(clc_upper);

%% 6) Plots
t = (0:N-1)/fs;
% a) Three leads (I, II, V1-like = leads 1,2,3)
f1 = figure('Color','w'); 
subplot(3,1,1); plot(t, ecg_noisy(1,:)); title('Lead 1'); xlim([0 6]);
subplot(3,1,2); plot(t, ecg_noisy(2,:)); title('Lead 2'); xlim([0 6]);
subplot(3,1,3); plot(t, ecg_noisy(3,:)); title('Lead 3'); xlim([0 6]);
xlabel('Time (s)'); 
saveas(f1, fullfile(outdir,"fig_ecg_three_leads.png"));

% b) Maternal vs fetal VCG loops (first 4 s)
idx = t<=4;
f2 = figure('Color','w'); 
plot3(vcg_m.x(idx), vcg_m.y(idx), vcg_m.z(idx)); hold on;
plot3(vcg_f.x(idx), vcg_f.y(idx), vcg_f.z(idx));
grid on; legend('Maternal VCG','Fetal VCG'); xlabel('X'); ylabel('Y'); zlabel('Z');
title('VCG Loops (first 4 s)');
saveas(f2, fullfile(outdir,"fig_vcg_loops.png"));

% c) RR histograms
f3 = figure('Color','w'); 
edges = 300:10:1200;
histogram(rr_m_ms, edges); hold on;
histogram(rr_f_ms, 150:5:600);
xlabel('RR (ms)'); ylabel('Count'); legend('Maternal','Fetal');
title('RR Interval Distributions');
saveas(f3, fullfile(outdir,"fig_rr_hist.png"));

%% 7) Console summary
fprintf('\n=== RESULTS ===\n');
fprintf('Mean SNR across 12 leads: %.1f dB\n', mean_snr_db);
fprintf('Beat-to-beat morphology variability (cosine dist, maternal-dom lead): %.3f\n', morph_var);
fprintf('HRV SDNN (ms): maternal %.1f, fetal %.1f\n', sdnn_m, sdnn_f);
fprintf('Avg cross-lead correlation: %.3f\n', crosslead_corr);
fprintf('Maternal beats detected: %d | Fetal beats detected: %d\n', numel(rpos_m), numel(rpos_f));

%% 8) Write Markdown report
rep = fullfile(outdir,"report.md");
fid = fopen(rep, 'w');
fprintf(fid, '# Multichannel ECG Modeling with Stochastic Extensions (OSET)\n\n');
fprintf(fid, '## Methods\n');
fprintf(fid, '- **Simulators:** vcg\\_gen\\_stochastic for maternal (%.0f bpm) and fetal (%.0f bpm) VCGs with per-beat jitter (HR, amplitude, width, phase).\n', f_m_bpm, f_f_bpm);
fprintf(fid, '- **Projection:** 3D VCG → 12-lead via a fixed 12×3 lead field (orthonormal). Replace with inverse-Dower for higher physiological realism.\n');
fprintf(fid, '- **Noise:** Baseline wander (0.2–0.5 Hz), EMG-like 20–50 Hz, and white noise to target %.0f dB SNR.\n', target_snr_db);
fprintf(fid, '- **Sampling:** fs=%d Hz, duration=%d s.\n\n', fs, dur_s);

fprintf(fid, '## Results\n');
fprintf(fid, '| Metric | Value |\n|---|---:|\n');
fprintf(fid, '| Mean SNR across 12 leads (dB) | %.1f |\n', mean_snr_db);
fprintf(fid, '| Morphology variability (cosine, maternal-dom lead) | %.3f |\n', morph_var);
fprintf(fid, '| HRV SDNN (ms) – maternal | %.1f |\n', sdnn_m);
fprintf(fid, '| HRV SDNN (ms) – fetal | %.1f |\n', sdnn_f);
fprintf(fid, '| Avg cross-lead correlation | %.3f |\n\n', crosslead_corr);

fprintf(fid, '### Figures\n');
fprintf(fid, '- ![](fig_ecg_three_leads.png)\n');
fprintf(fid, '- ![](fig_vcg_loops.png)\n');
fprintf(fid, '- ![](fig_rr_hist.png)\n\n');

fprintf(fid, '## Discussion: Why Stochastic Modeling Matters\n');
fprintf(fid, '- **Physiological realism:** Real ECGs exhibit micro-variability in RR intervals and morphology due to autonomic tone, respiration, and conduction variability. Adding controlled randomness avoids overly periodic "toy" signals.\n');
fprintf(fid, '- **Algorithm stress-testing:** Jittered morphology + noise exposes edge cases for QRS detection, segmentation, features, and classifiers, improving robustness and generalization.\n');
fprintf(fid, '- **Maternal–fetal use-case:** Two asynchronous stochastic sources create realistic overlap useful for fetal QRS detection and source separation benchmarking.\n');
fprintf(fid, '- **Extensibility:** The same framework supports alternans, state switching (PVCs), and time-warping to model broader pathologies and conditions.\n\n');

fprintf(fid, '## Limitations\n');
fprintf(fid, '- Lead field here is a proxy; use measured/inverse-Dower matrices for clinical realism.\n');
fprintf(fid, '- Fetal HRV is approximate with a crude high-pass detector; better separation/detectors will improve accuracy.\n');
fprintf(fid, '- Excessive parameter jitter can yield non-physiological shapes; tune within realistic ranges.\n');
fclose(fid);

fprintf('\nReport written to: %s\n', rep);
disp('Done.');

%% --------- Helper FUNCTIONS (kept inside file for portability) ---------
function ecg = project_vcg_to_ecg(vcg, L)
    % vcg: struct with x,y,z (1xN). L: Mx3.
    XYZ = [vcg.x(:) vcg.y(:) vcg.z(:)].';
    ecg = L * XYZ; % M x N
end

function y = add_ecg_noise(x, fs, snr_db, bw_amp, emg_scale)
    % x: MxN
    [M,N] = size(x);
    t = (0:N-1)/fs;

    % Baseline wander (sum of two very low sinusoids)
    bw = bw_amp*(0.7*sin(2*pi*0.33*t) + 0.3*sin(2*pi*0.22*t));
    bw = repmat(bw, M, 1);

    % EMG-like (band-limited 20–50 Hz)
    emg = randn(M,N);
    [b,a] = butter(2, [20 50]/(fs/2));
    for m=1:M
        emg(m,:) = filter(b,a,emg(m,:));
    end
    emg = emg * emg_scale;

    % White noise to hit target SNR
    sigp = mean(var(x,0,2));      % average signal power across leads
    np   = sigp / (10^(snr_db/10));
    wn   = sqrt(np) * randn(M,N);

    y = x + bw + emg + wn;
end

function rpos = simple_rpeaks(x, fs, varargin)
    % Crude R detector: bandpass + adaptive threshold + refractory
    p = inputParser; addParameter(p,'expected_bpm',[50 200]);
    parse(p, varargin{:});
    bpm = p.Results.expected_bpm;

    % bandpass 5–18 Hz (maternal-focused)
    xb = bandpass(x, [5 18], fs);

    % energy envelope
    env = movmean(abs(xb), round(0.050*fs));

    % threshold based on percentile
    thr = prctile(env, 92);
    cand = find(env > thr);

    % group by refractory (~0.25 s for maternal; adjust from bpm range)
    refr = round(0.25*fs);
    if numel(bpm)==2
        refr = round(60/max(bpm)*fs*0.6); % adapt to upper bpm
    end

    rpos = [];
    k = 1;
    while k <= numel(cand)
        win = cand(k):min(numel(env), cand(k)+refr-1);
        [~,ix] = max(env(win));
        rpos(end+1) = win(1)+ix-1; %#ok<AGROW>
        k = k + nnz(cand>=win(1) & cand<=win(end));
    end

    % sanity: remove peaks too close
    rpos = unique(rpos);
    rpos = rpos(:).';
end

function mv = beat_to_beat_variability(x, rpos, fs)
    % Compare consecutive beats by cosine distance on normalized beats
    if numel(rpos) < 3, mv = NaN; return; end
    pre = round(0.20*fs); post = round(0.10*fs);
    B = [];
    for i = 2:numel(rpos)-1
        s = max(1, rpos(i)-pre); e = min(length(x), rpos(i)+post);
        b = x(s:e);
        B = padcat_row(B, b);
    end
    if isempty(B), mv=NaN; return; end
    % z-normalize beats
    B = (B - mean(B,2))./max(std(B,0,2),1e-6);
    d = zeros(size(B,1)-1,1);
    for i=1:size(B,1)-1
        a = B(i,:); c = B(i+1,:);
        d(i) = 1 - (a*c.')/(norm(a)*norm(c)+eps);
    end
    mv = mean(d,'omitnan');
end

function M = padcat_row(M, row)
    if isempty(M), M = row; return; end
    len = size(M,2);
    if numel(row) == len
        M = [M; row];
    elseif numel(row) < len
        M = [M; [row, nan(1, len-numel(row))]];
    else
        M = [ [M, nan(size(M,1), numel(row)-len)]; row ];
    end
end

function R = axang2rotm(axang)
    % axang = [ux uy uz angle]; Rodrigues' rotation
    u = axang(1:3); th = axang(4);
    u = u/norm(u+eps);
    K = [  0   -u(3)  u(2);
          u(3)   0   -u(1);
         -u(2) u(1)    0 ];
    R = eye(3) + sin(th)*K + (1-cos(th))*(K*K);
end

end % main function

function [snr_per_lead_db, mean_snr_db] = snr_per_lead(clean_sig, noisy_sig)
% Compute SNR (dB) for each lead and the mean across leads
% clean_sig, noisy_sig: MxN matrices (M = leads, N = samples)
noise = noisy_sig - clean_sig;
sig_power = var(clean_sig, 0, 2);
noise_power = var(noise, 0, 2);
snr_per_lead_db = 10 * log10(sig_power ./ noise_power);
mean_snr_db = mean(snr_per_lead_db, 'omitnan');
end

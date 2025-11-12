
stochastic_ecg_metrics_demo()
function stochastic_ecg_metrics_demo()
% Main stochastic model + noise + metrics (SNR, HRV, morphology variability)

%% ---------------- Parameters ----------------
rng(123);
fs   = 500;                 % Hz
dur  = 30;                  % s
N    = fs*dur;
fHz  = 70/60;               % 70 BPM in Hz
theta0 = -2.1;
% Gaussian morphology (example)
alpha = struct('x',[0.85 0.50 0.30], 'y',[0.65 0.40 0.20], 'z',[0.75 0.50 0.24]);
b     = struct('x',[0.06 0.05 0.05], 'y',[0.06 0.05 0.05], 'z',[0.06 0.05 0.05]);
theta = struct('x',[-1.10 -0.15 0.22], 'y',[-1.05 -0.20 0.24], 'z',[-1.00 -0.25 0.26]);

% Stochastic deltas (beat-wise jitter)
f_deviations = 0.06;      % ±6% HR jitter
d_alpha = 0.10; d_b = 0.10; d_theta = 0.06;

% Projection to multi-lead (random orthonormal 12x3; swap with inverse-Dower if you have it)
A = randn(12,3); [Q,~] = qr(A,0); L = Q;

% Noise strengths
target_snr_db = 15;     % overall SNR target
bw_amp   = 0.05;        % baseline wander scale (~5% FS)
emg_gain = 0.02;        % EMG-like noise scale

%% ---------------- Synthesis ----------------
[vcg, phi] = vcg_gen_stochastic(N, fs, fHz, f_deviations, alpha, d_alpha, b, d_b, theta, d_theta, theta0); %#ok<ASGLU>
ecg_clean  = project_vcg_to_ecg(vcg, L);        % 12 x N

% Build noise and scale to target SNR
[ecg_noisy, parts] = add_noise_threeway(ecg_clean, fs, target_snr_db, bw_amp, emg_gain);

%% ---------------- Metrics ----------------
% SNR per lead
[snr_lead_db, snr_mean_db] = snr_per_lead(ecg_clean, ecg_noisy);

% Pick maternal-dominant lead by variance
[~, lead_idx] = max(var(ecg_clean,0,2));
x = ecg_noisy(lead_idx,:);

% R-peaks (simple but robust): bandpass + envelope + adaptive threshold
rpos = simple_rpeaks(x, fs);

% HRV: RR intervals → BPM mean & SDNN
rr = diff(rpos)/fs;                 % seconds
hr_bpm_inst = 60./rr;
hr_mean = mean(hr_bpm_inst);
sdnn_ms = std(rr)*1000;

% Morphology variability: cosine distance between consecutive normalized beats
morph_var = beat_to_beat_cosine_var(x, rpos, fs);

%% ---------------- Print & quick plots ----------------
fprintf('\n=== Stochastic ECG + Noise Metrics ===\n');
fprintf('Mean SNR across 12 leads: %.1f dB\n', snr_mean_db);
fprintf('HR mean: %.1f bpm | SDNN: %.1f ms | Morphology variability (cosine): %.3f\n', hr_mean, sdnn_ms, morph_var);

t = (0:N-1)/fs;
figure('Color','w'); 
subplot(3,1,1); plot(t, ecg_clean(lead_idx,:)); title(sprintf('Clean ECG (lead %d)', lead_idx)); xlim([0 6]); ylabel('mV');
subplot(3,1,2); plot(t, ecg_noisy(lead_idx,:)); hold on; stem(rpos/fs, 0.9*max(ecg_noisy(lead_idx,1:fs*6)) * ones(size(rpos)), 'r.'); xlim([0 6]);
title('Noisy ECG + detected R-peaks'); ylabel('mV');
subplot(3,1,3); plot(t, parts.baseline(lead_idx,:), t, parts.emg(lead_idx,:), t, parts.white(lead_idx,:));
xlim([0 6]); legend('baseline','EMG-like','white'); title('Noise components (first 6 s)'); xlabel('Time (s)');

end

%% ================== Helpers ==================

function ecg = project_vcg_to_ecg(vcg, L)
    XYZ = [vcg.x(:) vcg.y(:) vcg.z(:)].';
    ecg = L * XYZ; % M x N
end

function [y, parts] = add_noise_threeway(x, fs, target_snr_db, bw_amp, emg_gain)
% baseline wander (0.2–0.5 Hz), EMG-like (20–50 Hz), and white noise
    [M,N] = size(x); t = (0:N-1)/fs;
    parts = struct();

    % Baseline wander
    bw = bw_amp * (0.7*sin(2*pi*0.33*t) + 0.3*sin(2*pi*0.22*t));
    parts.baseline = repmat(bw, M,1);

    % EMG-like: band-limited 20–50 Hz (2nd-order Butterworth)
    emg = randn(M,N);
    [b1,a1] = butter(2, [20 50]/(fs/2));
    for m=1:M, emg(m,:) = filtfilt(b1,a1,emg(m,:)); end
    parts.emg = emg_gain * emg;

    % White noise (to be scaled)
    white = randn(M,N);
    parts.white = white; % temporary

    % Combine raw noises first (baseline + EMG), then scale white to hit target SNR
    n_partial = parts.baseline + parts.emg;
    sigp = mean(var(x,0,2));
    np_partial = mean(var(n_partial,0,2));

    % Determine white noise power needed
    desired_np = sigp / (10^(target_snr_db/10));
    np_white = max(desired_np - np_partial, 0);
    scale_white = sqrt(np_white);
    parts.white = scale_white .* parts.white;

    y = x + n_partial + parts.white;
end

function [snr_per_lead_db, mean_snr_db] = snr_per_lead(clean_sig, noisy_sig)
    noise = noisy_sig - clean_sig;
    sigp = var(clean_sig,0,2);
    np   = var(noise,0,2);
    snr_per_lead_db = 10*log10(sigp./np);
    mean_snr_db = mean(snr_per_lead_db,'omitnan');
end

function rpos = simple_rpeaks(x, fs)
% bandpass 5–18 Hz + envelope + adaptive threshold + refractory
    [b,a] = butter(2, [5 18]/(fs/2)); xb = filtfilt(b,a,x);
    env = movmean(abs(xb), round(0.050*fs));
    thr = prctile(env,92);
    cand = find(env > thr);

    refr = round(0.25*fs);     % ~250 ms refractory
    rpos = [];
    k = 1;
    while k <= numel(cand)
        win = cand(k):min(numel(env), cand(k)+refr-1);
        [~,ix] = max(env(win));
        rpos(end+1) = win(1)+ix-1; %#ok<AGROW>
        k = k + nnz(cand>=win(1) & cand<=win(end));
    end
    rpos = unique(rpos); rpos = rpos(:).';
    % basic sanity: keep peaks inside [0.2, dur-0.1] s
    rpos = rpos(rpos>round(0.2*fs) & rpos < numel(x)-round(0.1*fs));
end

function mv = beat_to_beat_cosine_var(x, rpos, fs)
% cosine distance between consecutive z-normalized beats
    if numel(rpos) < 3, mv = NaN; return; end
    pre = round(0.20*fs); post = round(0.10*fs);
    % collect equal-length beat snippets
    segs = cell(0,1);
    for i=2:numel(rpos)-1
        s = max(1, rpos(i)-pre); e = min(length(x), rpos(i)+post);
        segs{end+1} = x(s:e); %#ok<AGROW>
    end
    L = min(cellfun(@numel,segs));
    B = cell2mat(cellfun(@(v)v(1:L),segs,'uni',false)); % rows = beats
    % z-norm each beat
    B = (B - mean(B,2))./max(std(B,0,2),1e-6);
    % consecutive cosine distance
    d = zeros(size(B,1)-1,1);
    for k=1:numel(d)
        a = B(k,:); c = B(k+1,:);
        d(k) = 1 - (a*c.')/(norm(a)*norm(c)+eps);
    end
    mv = mean(d,'omitnan');
end

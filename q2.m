clear; clc; close all;

%% Parameters from McSharry & Clifford (TBME 2003) / ECGSYN
%% Part (a): Generate synthetic ECG (McSharry–Clifford model) and plot
clear; clc;
%
% ---------- Parameters (paper-aligned) ----------
mode_w   = 1;        % 1 = classic LF/HF RR process (paper); 2 = sinusoidal HRV (extension)
sfecg    = 256;      % output sampling rate (Hz)
N        = 256;      % ~number of beats
Anoise   = 0;     % mV noise
hrmean   = 60;       % bpm
hrstd    = 1;        % bpm (STD of HR)  -> drives RR variance
lfhfratio = 0.5;     % LF/HF power ratio (paper examples)
sfint    = 512;      % internal integration rate (Hz)

% PQRST morphology (paper Table I defaults)
ti_deg = [-70 -15 0 15 100];       % degrees, P Q R S T angles (relative to R)
ai     = [1.2 -5 30 -7.5 0.75];    % z "heights" (+/-)
bi     = [0.25 0.1 0.1 0.1 0.4];   % Gaussian widths

% ---------- Generate ECG using your ecgsyn ----------
[s, ipeaks, out] = ecgsyn(mode_w, sfecg, N, Anoise, hrmean, hrstd, lfhfratio, sfint, ti_deg, ai, bi);

% ---------- Time-domain plot ----------
t = (0:numel(s)-1)/sfecg;
figure('Color','w'); 
plot(t, s, 'LineWidth', 2,'color','k'); grid on;
xlabel('Time (s)'); ylabel('ECG (mV)');
title('Synthetic ECG (McSharry–Clifford model');

% Mark R-peaks (optional)
R = find(ipeaks==3);
hold on; plot(t(R), s(R), 'ro', 'MarkerSize', 4, 'DisplayName','R');
legend('ECG','R-peaks');
xlim([0.5,4.5])
saveas(gcf, 'synthetic_ecg.png');

%
%  Morphology phase-portrait: x vs z
z = out.X(:,3);     
figure('Color','w');
plot(out.X(:,1), z, 'k'); grid on;
xlabel('x'); ylabel('ECG z (mV)');
title('Morphology phase-portrait (x vs z)');
saveas(gcf, 'x_vs_z.png');
% Morphology phase-portrait: y vs z    
figure('Color','w');
plot(out.X(:,2), z, 'k'); grid on;
xlabel('y'); ylabel('ECG z (mV)');
title('Morphology phase-portrait (y vs z)');
saveas(gcf, 'y_vs_z.png');

% State-plane portrait: z vs dz/dt
dt = mean(diff(t));
dz = gradient(z, dt);
figure('Color','w');
plot(z, dz, 'k'); grid on;
xlabel('z (ECG)'); ylabel('dz/dt');
title('State-plane portrait (z vs dz/dt)');
saveas(gcf, 'z_vs_dz.png');

%

% --- Generate data (use your ecgsyn + derivsecgsyn) ---
[s, ipeaks, out] = ecgsyn(1, 256, 9, 0, 60, 1, 0.5, 512);  

% Use the raw z from the state (clean), then scale to match the paper look (~[-0.5, 1.5])
z  = out.X(:,3);
z  = 2 * (z - min(z)) / (max(z) - min(z)) - 0.5;   % [-0.5, 1.5] dynamic range

% --- 3D trajectory ---
figure; hold on;
plot3(out.X(:,1), out.X(:,2), z, 'k', 'LineWidth', 2);   % main thick curve

% --- dashed baseline orbit (z ~ constant) ---
th  = linspace(-pi, pi, 400);
z0  = median(z);                                         % baseline level
plot3(cos(th), sin(th), z0*ones(size(th)), 'k--', 'LineWidth', 1.0);

% --- optional dotted markers on the baseline circle (like the figure) ---
thdots = linspace(out.ti(3)-0.3, out.ti(3)+0.3, 3);      % around the R region
plot3(cos(thdots), sin(thdots), z0*ones(size(thdots)), 'k.', 'MarkerSize', 18);

% --- P, Q, R, S, T landmarks ---
labs = {'P','Q','R','S','T'};
for i = 1:numel(out.ti)
    % nominal (x,y) at angle ti(i)
    xi = cos(out.ti(i));
    yi = sin(out.ti(i));
    % find nearest actual sample to get the corresponding z value
    [~, idx] = min(abs(wrapToPi(out.theta - out.ti(i))));
    zi = z(idx);
    plot3(xi, yi, zi, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
    text(xi, yi, zi + 0.06, labs{i}, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end


xlabel('x'); ylabel('y'); zlabel('z');
title('3D ECG Trajectory with PQRST Landmarks');
axis equal; grid on; box on;
xlim([-1 1]); ylim([-1 1]); zlim([-0.5 1.6]);

% choose a view 
view(35, 20);    
saveas(gcf, 'xyz.png');


%% Part B 
% Parameter Variation Analysis for ALL P,Q,R,S,T (ai & bi) — FULL ECG PLOTS
% Uses your ecgsyn() and derivsecgsyn() functions (already defined).
% Plots complete ECG for each variation, not just one beat.

clear; clc;

% ---------------- Simulation Settings ----------------
mode_w   = 2;         % sinusoidal mode, but depth=0 => constant HR
sfecg    = 256;       % ECG sampling frequency [Hz]
N        = 300;       % number of beats (long ECG)
Anoise   = 0.00;      % no added noise
hrmean   = 60;        % bpm
depth    = 0.00;      % no HRV
fmod     = 0.10;      % Hz (irrelevant when depth=0)
sfint    = 512;       % internal ODE sampling freq [Hz]

% PQRST parameters from paper
ti_deg = [-70 -15 0 15 100];          % angular positions
ai0    = [1.2 -5 30 -7.5 0.75];       % amplitudes (P Q R S T)
bi0    = [0.25 0.1 0.1 0.1 0.4];      % widths (P Q R S T)
labels = {'P','Q','R','S','T'};

% Variation scales
scaleA = [0.6 0.8 1.0 1.2 1.5];
scaleB = [0.6 0.8 1.0 1.2 1.5];

% Helper handles
gen = @(ai,bi) ecgsyn(mode_w, sfecg, N, Anoise, hrmean, depth, fmod, sfint, ti_deg, ai, bi);
t_of = @(s) (0:numel(s)-1)/sfecg;

% Output folder for saving figures
outdir = 'ECG_Variations';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

% ---------- Reference ECG ----------
[s_ref, ~] = gen(ai0, bi0);
t = t_of(s_ref);

fig = figure('Color','w');
plot(t, s_ref, 'k', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('ECG (mV)');
title(sprintf('Reference ECG — %d beats at %g bpm', N, hrmean));
grid on;
xlim([0.5 4.5]);
saveas(fig, fullfile(outdir, 'Reference_ECG.png'));

% ---------- VARIATION OF ai (AMPLITUDE) ----------
for k = 1:5
    fig = figure('Color','w'); hold on; grid on;
    for sc = 1:numel(scaleA)
        ai = ai0; ai(k) = ai0(k) * scaleA(sc);
        [s, ~] = gen(ai, bi0);
        t = t_of(s);
        plot(t, s, 'LineWidth', 1, ...
             'DisplayName', sprintf('a_{%s} × %.1f', labels{k}, scaleA(sc)));
    end
    xlabel('Time (s)');
    ylabel('ECG (mV)');
    title(sprintf('Effect of a_{%s} (Amplitude) on Full ECG', labels{k}));
    legend('show','Location','best');
    xlim([0.5 4.5]);
    fname = sprintf('Effect_of_a_%s.png', labels{k});
    saveas(fig, fullfile(outdir, fname));
end

% ---------- VARIATION OF bi (WIDTH) ----------
for k = 1:5
    fig = figure('Color','w'); hold on; grid on;
    for sc = 1:numel(scaleB)
        bi = bi0; bi(k) = bi0(k) * scaleB(sc);
        [s, ~] = gen(ai0, bi);
        t = t_of(s);
        plot(t, s, 'LineWidth', 1, ...
             'DisplayName', sprintf('b_{%s} × %.1f', labels{k}, scaleB(sc)));
    end
    xlabel('Time (s)');
    ylabel('ECG (mV)');
    title(sprintf('Effect of b_{%s} (Width) on Full ECG', labels{k}));
    legend('show','Location','best');
    xlim([0.5 4.5]);
    fname = sprintf('Effect_of_b_%s.png', labels{k});
    saveas(fig, fullfile(outdir, fname));
end

fprintf('All figures saved in folder: "%s"\n', outdir);

%% Part C 
% === Parameters ===
mode_w   = 2;         % 1 = classic RR-driven, 2 = sinusoidal HRV (your iii)
sfecg    = 256;       % Hz (output sampling rate)
N        = 256;       % approx. number of beats
Anoise   = 0;      % mV of uniform noise
hrmean   = 60;        % bpm
depth    = 0.40;      % ±10% HR modulation (this maps to hrstd in mode 2)
fmod     = 0.10;      % Hz HRV modulation (this maps to lfhfratio in mode 2)
sfint    = 512;       % Hz (internal)
ti_deg   = [-70 -15 0 15 100];              % default morphology
ai       = [1.2 -5 30 -7.5 0.75];
bi       = [0.25 0.1 0.1 0.1 0.4];

% === Generate ECG with sinusoidal HRV ===
% In mode 2: hrstd := depth (fraction), lfhfratio := fmod (Hz)
[s, ipeaks] = ecgsyn(mode_w, sfecg, N, Anoise, hrmean, depth, fmod, sfint, ti_deg, ai, bi);

% === Time vector and basic plot ===
t = (0:numel(s)-1)/sfecg;
figure; plot(t, s, 'LineWidth', 1);
xlabel('Time (s)'); ylabel('ECG (mV)');
title(sprintf('Synthetic ECG with Sinusoidal HRV: depth=%.0f%%, f_{mod}=%.2f Hz', 100*depth, fmod));
grid on;
xlim([0 20])
saveas(gcf, 'HRV.png');
% % === Derive R-peaks, RR intervals, and instantaneous HR (for visualization) ===
% R = find(ipeaks==3);
% R = R(:);
% RR_sec = diff(R)/sfecg;           % seconds per beat
% HR_bpm = 60 ./ RR_sec;

% % Plot instantaneous HR vs time (located between beats)
% tHR = t(R(1:end-1)) + RR_sec/2;
% figure; plot(tHR, HR_bpm, 'o-','LineWidth',1);
% xlabel('Time (s)'); ylabel('Instantaneous HR (bpm)');
% title('Instantaneous Heart Rate from R–R intervals');
% grid on;



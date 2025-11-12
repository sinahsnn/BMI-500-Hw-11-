function demo_vcg_generators()
% Demo for OSET VCG generators (multichannel, 3D) → 12-lead ECG projection.
% Runs:
%   1) vcg_gen_abnormal
%   2) vcg_gen_direct_sum
%   3) vcg_gen_state_space
%   4) vcg_gen_stochastic
%   5) vcg_gen_var_hr  (variable HR, length may differ from others)
%
% It saves simple figures for each and prints a few basics.

%% Common parameters
rng(42);
fs   = 500;             % Hz
dur  = 10;              % seconds (for fixed-length generators)
N    = dur*fs;          % samples
fHz  = 70/60;           % 70 bpm in Hz
theta0 = -2.1;

% Gaussian parameters (same base morphology across demos)
alpha = struct('x',[0.85 0.50 0.30], 'y',[0.65 0.40 0.20], 'z',[0.75 0.50 0.24]);
b     = struct('x',[0.06 0.05 0.05], 'y',[0.06 0.05 0.05], 'z',[0.06 0.05 0.05]);
theta = struct('x',[-1.10 -0.15 0.22], 'y',[-1.05 -0.20 0.24], 'z',[-1.00 -0.25 0.26]);

% Lead field for 12-lead projection (random orthonormal 12x3; swap with inverse-Dower if you have it)
A = randn(12,3); [Q,~] = qr(A,0); L = Q;

% Output directory
outdir = "vcg_demo_outputs"; if ~exist(outdir,"dir"), mkdir(outdir); end

%% 1) vcg_gen_abnormal — two-state alternation (e.g., T-wave alternans style)
fprintf('\n[1] vcg_gen_abnormal...\n');
alphaA = alpha;
alphaB = alpha; alphaB.x(3)=alphaB.x(3)*0.75; alphaB.y(3)=alphaB.y(3)*0.75; alphaB.z(3)=alphaB.z(3)*0.75; % smaller T
alpha_states = [alphaA alphaB];
b_states     = [b b];
theta_states = [theta theta];
STM = [0 1; 1 0];   % strict alternation A<->B
S0  = 1;

[vcg_abn, phi_abn] = vcg_gen_abnormal(N, fs, fHz, alpha_states, b_states, theta_states, theta0, STM, S0); %#ok<ASGLU>
ecg_abn = project_vcg_to_ecg(vcg_abn, L);
quick_plot(ecg_abn, fs, 6, 'vcg\_gen\_abnormal → 12-lead (first 6 s)');
saveas(gcf, fullfile(outdir, "fig_abnormal.png"));

%% 2) vcg_gen_direct_sum — vectorized direct sum
fprintf('[2] vcg_gen_direct_sum...\n');
[vcg_dir, phi_dir] = vcg_gen_direct_sum(N, fs, fHz, alpha, b, theta, theta0); %#ok<ASGLU>
ecg_dir = project_vcg_to_ecg(vcg_dir, L);
quick_plot(ecg_dir, fs, 6, 'vcg\_gen\_direct\_sum → 12-lead (first 6 s)');
saveas(gcf, fullfile(outdir, "fig_direct_sum.png"));

%% 3) vcg_gen_state_space — differential/state-space update
fprintf('[3] vcg_gen_state_space...\n');
[vcg_ss, phi_ss] = vcg_gen_state_space(N, fs, fHz, alpha, b, theta, theta0); %#ok<ASGLU>
ecg_ss = project_vcg_to_ecg(vcg_ss, L);
quick_plot(ecg_ss, fs, 6, 'vcg\_gen\_state\_space → 12-lead (first 6 s)');
saveas(gcf, fullfile(outdir, "fig_state_space.png"));

%% 4) vcg_gen_stochastic — beat-wise random HR & morphology
fprintf('[4] vcg_gen_stochastic...\n');
f_deviations = 0.06;     % ±6% beat-to-beat HR jitter
d_alpha = 0.10; d_b = 0.10; d_theta = 0.06;
[vcg_sto, phi_sto] = vcg_gen_stochastic(N, fs, fHz, f_deviations, alpha, d_alpha, b, d_b, theta, d_theta, theta0); %#ok<ASGLU>
ecg_sto = project_vcg_to_ecg(vcg_sto, L);
quick_plot(ecg_sto, fs, 6, 'vcg\_gen\_stochastic → 12-lead (first 6 s)');
saveas(gcf, fullfile(outdir, "fig_stochastic.png"));

%% 5) vcg_gen_var_hr — time-varying HR trajectory (MATRIX)
fprintf('[5] vcg_gen_var_hr (MATRIX)...\n');
% Build a per-beat HR vector around 150 bpm with gentle modulation
meanHR = 150; modAmp = 8; numBeats = ceil(dur * meanHR/60);   % approx beats in 'dur' seconds
HR_vec = meanHR + modAmp*sin(2*pi*(1:numBeats)/10);           % BPM per beat (slow sinus modulation)
theta0_f = 0.8;

[vcg_vhr, phi_vhr] = vcg_gen_var_hr(HR_vec, fs, alpha, b, theta, theta0_f, 'MATRIX'); %#ok<ASGLU>
ecg_vhr = project_vcg_to_ecg(vcg_vhr, L);
fprintf('   var_hr output length: %d samples (may differ from N=%d)\n', size(ecg_vhr,2), N);

% Plot first 6 s or full if shorter
maxs = min(6, size(ecg_vhr,2)/fs);
quick_plot(ecg_vhr, fs, maxs, 'vcg\_gen\_var\_hr → 12-lead (first seconds)');
saveas(gcf, fullfile(outdir, "fig_var_hr.png"));

%% Console notes
fprintf('\nSaved figures in folder: %s\n', outdir);
fprintf('Done.\n');

%% ----------------- helpers -----------------
function ecg = project_vcg_to_ecg(vcg, L)
    % vcg: struct with fields x,y,z (1xN). L: Mx3 lead field.
    XYZ = [vcg.x(:) vcg.y(:) vcg.z(:)].';
    ecg = L * XYZ; % M x N
end

function quick_plot(ecgM, fs, secs, ttl)
    t = (0:size(ecgM,2)-1)/fs;
    show = t <= secs;
    figure('Color','w'); 
    subplot(3,1,1); plot(t(show), ecgM(1,show)); title([ttl ' — Lead 1']); ylabel('mV'); grid on;
    subplot(3,1,2); plot(t(show), ecgM(2,show)); title('Lead 2'); ylabel('mV'); grid on;
    subplot(3,1,3); plot(t(show), ecgM(3,show)); title('Lead 3'); ylabel('mV'); xlabel('Time (s)'); grid on;
end

end

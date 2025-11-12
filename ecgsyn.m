function [s, ipeaks, out] = ecgsyn(mode_w,sfecg,N,Anoise,hrmean,hrstd,lfhfratio,sfint,ti,ai,bi)
% [s, ipeaks] = ecgsyn(mode_w,sfecg,N,Anoise,hrmean,hrstd,lfhfratio,sfint,ti,ai,bi)
% Produces synthetic ECG:
%   s: ECG (mV)
%   ipeaks: labels for PQRST peaks: P(1), Q(2), R(3), S(4), T(5)
%
% mode_w:
%   1 -> classic model: omega(t) from RR(t) process (LF/HF HRV)
%   2 -> sinusoidal HRV: omega(t) from HR(t) = hrmean*(1 + depth*sin(2*pi*fmod*t))
%        (in this mode, hrstd := depth (fraction), lfhfratio := fmod (Hz))
%
% Default values (if omitted):
%   sfecg=256, N=256, Anoise=0, hrmean=60, hrstd=1 (or 0.10 depth in mode 2),
%   lfhfratio=0.5 (or 0.10 Hz fmod in mode 2), sfint=512,
%   ti=[-70 -15 0 15 100] deg, ai=[1.2 -5 30 -7.5 0.75], bi=[0.25 0.1 0.1 0.1 0.4]

% --- defaults (respecting new arg order) ---
if nargin < 1 || isempty(mode_w), mode_w = 1; end
if nargin < 2 || isempty(sfecg),  sfecg  = 256; end
if nargin < 3 || isempty(N),      N      = 256; end
if nargin < 4 || isempty(Anoise), Anoise = 0;   end
if nargin < 5 || isempty(hrmean), hrmean = 60;  end
if nargin < 6 || isempty(hrstd),  hrstd  = (mode_w==1)*1 + (mode_w~=1)*0.10; end % 1 bpm or 0.10 depth
if nargin < 7 || isempty(lfhfratio), lfhfratio = (mode_w==1)*0.5 + (mode_w~=1)*0.10; end % 0.5 or 0.10 Hz
if nargin < 8 || isempty(sfint),  sfint  = 512; end
if nargin < 9 || isempty(ti),     ti     = [-70 -15 0 15 100]; end
ti = ti*pi/180;  % radians
if nargin <10 || isempty(ai),     ai     = [1.2 -5 30 -7.5 0.75]; end
if nargin <11 || isempty(bi),     bi     = [0.25 0.1 0.1 0.1 0.4]; end

% --- adjust extrema for mean HR (paper's scaling) ---
hrfact  = sqrt(hrmean/60);
hrfact2 = sqrt(hrfact);
bi = hrfact*bi;
ti = [hrfact2 hrfact 1 hrfact hrfact2].*ti;

% --- check sfint multiple of sfecg ---
q  = round(sfint/sfecg);
qd = sfint/sfecg;
if q ~= qd
   error(['Internal sampling frequency (sfint) must be an integer multiple ' ... 
'of the ECG sampling frequency (sfecg). Your current choices are: ' ... 
'sfecg = ' int2str(sfecg) ' and sfint = ' int2str(sfint) '.']);
end

% --- print config ---
fid = 1;
fprintf(fid,'ECG sampled at %d Hz\n',sfecg);
fprintf(fid,'Approximate number of heart beats: %d\n',N);
fprintf(fid,'Measurement noise amplitude: %g mV\n',Anoise);
fprintf(fid,'mode_w = %d\n',mode_w);
fprintf(fid,'Heart rate mean: %g bpm\n',hrmean);
if mode_w==1
    fprintf(fid,'Heart rate std: %g bpm\n',hrstd);
    fprintf(fid,'LF/HF ratio: %g\n',lfhfratio);
else
    fprintf(fid,'Sinusoidal depth (fraction): %g\n',hrstd);
    fprintf(fid,'Sinusoidal fmod (Hz): %g\n',lfhfratio);
end
fprintf(fid,'Internal sampling frequency: %g\n',sfint);
fprintf(fid,'      P  Q  R  S  T\n'); 
fprintf(fid,'ti = [%g %g %g %g %g] radians\n',ti(1),ti(2),ti(3),ti(4),ti(5));
fprintf(fid,'ai = [%g %g %g %g %g]\n',ai(1),ai(2),ai(3),ai(4),ai(5));
fprintf(fid,'bi = [%g %g %g %g %g]\n',bi(1),bi(2),bi(3),bi(4),bi(5));

% --- build time base and RR / parameters depending on mode ---
dt = 1/sfint;

if mode_w==1
    % LF/HF HRV (original)
    flo = 0.1; fhi = 0.25; flostd = 0.01; fhistd = 0.01;

    % RR synthesis horizon (seconds): ~N beats at mean RR
    rrmean = 60/hrmean;
    sampfreqrr = 1;
    trr = 1/sampfreqrr; 
    Nrr = 2^(ceil(log2(N*rrmean/trr)));

    % RR(t) at 1 Hz, then upsample to sfint Hz
    rr0 = rrprocess(flo,fhi,flostd,fhistd,lfhfratio,hrmean,hrstd,sampfreqrr,Nrr);
    rr  = interp(rr0, sfint);

    % piecewise-constant rrn aligned to integration grid
    rrn = zeros(length(rr),1);
    tecg = 0; i = 1;
    while i <= length(rr)
        tecg = tecg + rr(i);
        ip = round(tecg/dt);
        rrn(i:ip) = rr(i);
        i = ip + 1;
    end
    Nt = ip;                                 % number of internal steps
    Tspan = 0:dt:(Nt-1)*dt;
    rr_arg = rrn;                            % pass RR array to ODE
else
    % Sinusoidal HRV (iii): use hrstd as depth, lfhfratio as fmod
    depth = hrstd;           % e.g., 0.10 for ±10%
    fmod  = lfhfratio;       % Hz, e.g., 0.10

    % choose duration to cover ~N beats at mean RR
    total_time = N * (60/hrmean);           % seconds
    Nt = ceil(total_time * sfint);
    Tspan = 0:dt:(Nt-1)*dt;

    % pack parameters for derivsecgsyn (mode 2 path)
    rr_arg = [hrmean, depth, fmod];
end

% --- integrate system using ODE45 ---
fprintf(fid,'Integrating dynamical system\n');
x0 = [1,0,0.04];
[T,X0] = ode45('derivsecgsyn', Tspan, x0, [], rr_arg, sfint, ti, ai, bi, mode_w);

% --- downsample to sfecg ---
Xds  = X0(1:q:end,:);          % downsampled states [x y z] @ sfecg
Tds  = Tspan(1:q:end).';       % downsampled time (column)

% keep legacy variable name for detectpeaks (expects X)
X    = Xds;

% --- phase-plane helper signals ---
Tint       = Tspan(:);         % internal time @ sfint
Xint       = X0;               % internal states @ sfint
theta_int  = atan2(Xint(:,2), Xint(:,1));
theta      = atan2(Xds(:,2),  Xds(:,1));

% --- package outputs (add out as 3rd output of ecgsyn) ---
out.Tint      = Tint;
out.Xint      = Xint;
out.theta_int = theta_int;

out.T         = Tds;
out.X         = Xds;
out.theta     = theta;

out.ti        = ti;                    % P,Q,R,S,T angles (radians)
out.labels    = {'P','Q','R','S','T'};




% --- detect peaks ---
ipeaks = detectpeaks(X, ti, sfecg);

% --- scale ECG to ~[-0.4, 1.2] mV and add noise ---
z = X(:,3);
zmin = min(z); zmax = max(z); zrange = zmax - zmin;
z = (z - zmin)*(1.6)/zrange - 0.4;
eta = 2*rand(length(z),1)-1;
s = z + Anoise*eta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rr = rrprocess(flo, fhi, flostd, fhistd, lfhfratio, hrmean, hrstd, sfrr, n)
w1 = 2*pi*flo;  w2 = 2*pi*fhi;
c1 = 2*pi*flostd; c2 = 2*pi*fhistd;
sig2 = 1; sig1 = lfhfratio;
rrmean = 60/hrmean;
rrstd  = 60*hrstd/(hrmean*hrmean);

df = sfrr/n;
w = (0:n-1)'*2*pi*df;
dw1 = w-w1; dw2 = w-w2;

Hw1 = sig1*exp(-0.5*(dw1/c1).^2)/sqrt(2*pi*c1^2);
Hw2 = sig2*exp(-0.5*(dw2/c2).^2)/sqrt(2*pi*c2^2);
Hw = Hw1 + Hw2;
Hw0 = [Hw(1:n/2); Hw(n/2:-1:1)];
Sw = (sfrr/2)*sqrt(Hw0);

ph0 = 2*pi*rand(n/2-1,1);
ph = [0; ph0; 0; -flipud(ph0)];
SwC = Sw .* exp(1j*ph);
x = (1/n)*real(ifft(SwC));

xstd = std(x);
ratio = rrstd/xstd;
rr = rrmean + x*ratio;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ind = detectpeaks(X, thetap, sfecg)
N = length(X);
theta = atan2(X(:,2),X(:,1));
ind0 = zeros(N,1);
for i=1:N-1
   a = ( (theta(i) <= thetap) & (thetap <= theta(i+1)) );
   j = find(a==1);
   if ~isempty(j)
      d1 = thetap(j) - theta(i);
      d2 = theta(i+1) - thetap(j);
      if d1 < d2, ind0(i) = j; else, ind0(i+1) = j; end
   end
end

d = ceil(sfecg/64);
d = max([2 d]);   % (left un-suppressed to match original behavior)
ind = zeros(N,1);
z = X(:,3);
zmin = min(z); zmax = max(z);
zext = [zmin zmax zmin zmax zmin];
sext = [1 -1 1 -1 1];
for i=1:5
   ind1 = find(ind0==i); n = length(ind1);
   Z = ones(n,2*d+1)*zext(i)*sext(i);
   for j=-d:d
      k = find( (1 <= ind1+j) & (ind1+j <= N) );
      Z(k,d+j+1) = z(ind1(k)+j)*sext(i);
   end
   [~, ivmax] = max(Z,[],2);
   iext = ind1 + ivmax-d-1;
   ind(iext) = i;
end

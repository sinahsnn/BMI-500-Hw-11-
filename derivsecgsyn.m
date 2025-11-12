function dxdt = derivsecgsyn(t,x,flag,rr,sfint,ti,ai,bi,mode_w)
% dxdt = derivsecgsyn(t,x,flag,rr,sfint,ti,ai,bi,mode_w)
% ODE file for generating the synthetic ECG
% rr: 
%   - if mode_w==1: RR(t) process sampled at sfint (seconds/beat)
%   - otherwise   : parameter vector [hrmean, depth, fmod]
% sfint: Internal sampling frequency [Hz]
% ti: angles of extrema [rad], ai: z-positions, bi: widths
%
% mode_w:
%   1 -> classic model, omega(t) = 2*pi / RR(t)
%   2 -> (or anything else) sinusoidal HRV: omega(t) from HR(t) = hrmean*(1+depth*sin(2*pi*fmod*t))

% --- geometry & radial dynamics (unchanged) ---
ta = atan2(x(2),x(1));
r0 = 1;
a0 = 1.0 - sqrt(x(1)^2 + x(2)^2)/r0;

% --- angular velocity selection ---
if mode_w==1
    % Use RR array (seconds/beat) at internal rate sfint
    ip = 1 + floor(t*sfint);
    % guard against bounds and nonpositive RR
    ip = max(1, min(ip, numel(rr)));
    rr_ip = max(rr(ip), 1e-3);
    w0 = 2*pi/rr_ip;  % rad/s
else
    % Sinusoidal HRV: rr argument carries [hrmean, depth, fmod]
    if ~isempty(rr)
        if numel(rr) >= 3
            hrmean = rr(1); depth = rr(2); fmod = rr(3);
        else
            hrmean = rr(1); depth = 0.10; fmod = 0.10;
        end
    else
        hrmean = 60; depth = 0.10; fmod = 0.10;
    end
    % instantaneous heart rate (bpm), clamped to avoid zero/negative
    hr_inst = max(hrmean * (1 + depth * sin(2*pi*fmod*t)), 5);
    % omega(t) = 2*pi * HR(t)/60  [rad/s]
    w0 = 2*pi*(hr_inst/60);
end

% --- baseline wander (respiratory-like) ---
fresp = 0.25;
zbase = 0.005*sin(2*pi*fresp*t);

% --- core ODEs ---
dx1dt = a0*x(1) - w0*x(2);
dx2dt = a0*x(2) + w0*x(1);

dti = rem(ta - ti, 2*pi);
dx3dt = - sum(ai.*dti.*exp(-0.5*(dti./bi).^2)) - 1.0*(x(3) - zbase);

dxdt = [dx1dt; dx2dt; dx3dt];
end

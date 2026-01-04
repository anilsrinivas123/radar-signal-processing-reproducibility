function est = candan_three_point_freq(x, fs)
% CANDAN_THREE_POINT_FREQ  Fine frequency estimate from 3 DFT samples
%   est = candan_three_point_freq(x,fs)
%   x : real or complex data vector (length N)
%   fs: sampling frequency in Hz
%   est is a struct with fields:
%     N       - length of DFT used
%     kp      - peak bin index (0-based)
%     delta   - fractional bin offset (units of bins, in [-0.5,0.5] typically)
%     omega   - estimated frequency (radians/sample)
%     freqHz  - estimated frequency in Hz
%
% Example:
%   fs = 8000;
%   N = 1024;
%   t = (0:N-1)'/fs;
%   x = exp(1j*(2*pi*1234.567*t));    % complex tone at 1234.567 Hz
%   est = candan_three_point_freq(x, fs);
%   fprintf('freq = %.6f Hz\n', est.freqHz);

if nargin<2, fs = 1; end

% ensure column
x = x(:);
N = length(x);

% if input is real, convert to complex analytic signal to avoid mirror peaks
if isreal(x)
    x = hilbert(x);
end

% compute N-point DFT (use FFT)
R = fft(x, N);

% magnitude and peak bin
[~, kmax] = max(abs(R));
% MATLAB indices: 1..N. Convert to 0-based for formula clarity
kp0 = kmax - 1;

% helper to wrap indices (0-based)
wrap = @(k) mod(k, N);

% extract R[kp-1], R[kp], R[kp+1] (in 0-based then convert to MATLAB 1-based)
R_m1 = R(wrap(kp0-1)+1);
R_0  = R(wrap(kp0)+1);
R_p1 = R(wrap(kp0+1)+1);

% compute correction factor gamma = tan(pi/N)/(pi/N)
gamma = tan(pi/N) / (pi/N);

% compute fraction (complex). Use denominator guard for numerical stability.
den = (2*R_0 - R_m1 - R_p1);
if abs(den) < eps
    frac = 0;
else
    frac = (R_m1 - R_p1) / den;
end

% Candan estimator uses the real part
delta_hat = gamma * real(frac);

% clamp to -0.5..0.5 (optional but practical)
if delta_hat > 0.5, delta_hat = 0.5; end
if delta_hat < -0.5, delta_hat = -0.5; end

% estimated frequency (radians/sample) and Hz
omega_hat = 2*pi*(kp0 + delta_hat)/N;        % radians/sample
freqHz    = omega_hat/(2*pi) * fs;         % Hz

% return results in struct
%est.N = N;
%est.kp = kp0;
%est.delta = delta_hat;
%est.omega = omega_hat;
%est.freqHz = freqHz;
est = freqHz;
%est.R = [R_m1, R_0, R_p1];  % (optional) include the three DFT samples

end

function [f_est, f_crlb]   = freq_est_fft_parabolic(y, fs, a, sigma2)
% Frequency estimation via FFT + parabolic interpolation
% Also computes CRLB for frequency estimation
%
% Inputs:
%   y      : input real signal (length N)
%   fs     : sampling frequency (Hz)
%   a      : signal amplitude (for CRLB)
%   sigma2 : noise variance (for CRLB)
%
% Outputs:
%   f_est  : estimated frequency (Hz)
%   f_crlb : CRLB on frequency estimate variance (Hz^2)

N = length(y);

% === 1) Apply window (Hamming for leakage reduction) ===
beta = 0.5;
w = kaiser(N, beta).';
%w = hamming(N).';
yw = y(:).' .* w;  % row vector

% === 2) Zero-pad and FFT ===
L = 2*N;   % zero-padding factor (adjustable)
Y = fft(yw, L);
[~, kmax] = max(abs(Y));  % coarse FFT bin index

% === 3) Parabolic interpolation ===
% Wrap-around indexing for neighbors
km1 = mod(kmax-2, L) + 1;
kp1 = mod(kmax,   L) + 1;

Pk_m1 = abs(Y(km1));
Pk    = abs(Y(kmax));
Pk_p1 = abs(Y(kp1));

delta = 0.5 * (Pk_m1 - Pk_p1) / (Pk_m1 - 2*Pk + Pk_p1);

k_refined = kmax - 1 + delta;  % fractional bin index (0-based)

% === 4) Frequency estimate ===
f_est = (k_refined / L) * fs;  % Hz

% Wrap to Nyquist interval
if f_est > fs/2
    f_est = f_est - fs;
end

% === 5) CRLB for frequency estimation ===
% Model: y[n] = a*cos(2π f n + phi) + w[n], w ~ N(0,σ²)
% CRLB(var(f_hat)) = 6*sigma² / (pi^2 * a^2 * N*(N^2 - 1))   [cycles/sample^2]
%
% Reference: Kay, "Fundamentals of Statistical Signal Processing: Estimation Theory"
%
% Convert to Hz^2: multiply by fs^2
f_crlb = (6 * sigma2 * fs^2) / (pi^2 * a^2 * N * (N^2 - 1));

end

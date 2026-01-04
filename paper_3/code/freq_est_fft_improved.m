function f_est = freq_est_fft_improved(y, fs)
% Frequency estimation via FFT + Jacobsen complex interpolation
% Computes correct CRLB for the used window

% Input Validation
N = length(y);
if N < 3
    error('Input signal y must have length >= 3.');
end

% === 1) Apply window (Kaiser) ===
beta = 0.5;
w = kaiser(N, beta);
w = w / norm(w); % Optional: normalize window energy
n = (0:N-1)';
yw = y(:) .* w;

% Calculate window efficiency factor (eta) for CRLB
% The term (n - (N-1)/2) centers the time vector to minimize numerical issues
%eta = ( (w.' * (n - (N-1)/2) ) ^2 ) / ( N * (w' * w) );

% === 2) Zero-pad and FFT ===
L = 2^nextpow2(4*N); % More aggressive padding for interpolation
Y = fft(yw, L);
[~, kmax] = max(abs(Y(1:floor(L/2)))); % Find peak in first half (real signal)

% === 3) Jacobsen (Complex) Interpolation ===
% Adjust indices for 1-based indexing. Y(kmax) is our peak.
% We use samples at kmax, and its neighbors kmax-1 and kmax+1.
% Handle boundaries:
k1 = max(kmax-1, 1);
k2 = kmax;
k3 = min(kmax+1, L);

Y1 = Y(k1);
Y2 = Y(k2);
Y3 = Y(k3);

% Jacobsen's formula
d = -real( (Y3 - Y1) / (2*Y2 - Y1 - Y3) );
k_refined = (k2 - 1) + d; % 'k2-1' converts to 0-index, then add delta

% === 4) Frequency estimate ===
f_est = (k_refined / L) * fs;

% Wrap to Nyquist interval (-fs/2, fs/2]
f_est = mod(f_est + fs/2, fs) - fs/2;
end
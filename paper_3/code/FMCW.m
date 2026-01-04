clc; clear; close all;

%% Parameters
c = 3e8;              % Speed of light (m/s)
fc = 77e9;            % Carrier frequency (Hz)
lambda = c / fc;      % Wavelength (m)

B = 200e6;            % Chirp bandwidth (Hz)
T = 20e-6;            % Chirp duration (s)
fs = 2 * B;           % Sampling frequency (Hz)
N = round(T * fs);    % Number of samples per chirp
t = (0:N-1)/fs;       % Time vector

%% Target Parameters
R_true = 100;         % Range to target (m)
v_true = 30;          % Radial velocity (m/s)

tau = 2 * R_true / c;                  % Time delay
fd = 2 * v_true / lambda;              % Doppler shift
doppler_scale = 1 - 2*v_true/c;        % Time scaling approx (for baseband)

%% Transmitted Chirp (baseband)
mu = B / T;  % Chirp rate
x = cos(pi * mu * t.^2);  % Baseband LFM chirp

%% Received Echo (delayed and Doppler scaled)
t_echo = (t - tau) / doppler_scale;    % Apply delay and Doppler time scaling
y = zeros(1, N);
valid_idx = t_echo > 0 & t_echo < T;
y(valid_idx) = cos(pi * mu * t_echo(valid_idx).^2);

%% Matched Filter (correlate received signal with transmitted)
z = xcorr(y, x);   % Cross-correlation
lags = (-N+1):(N-1);
[~, peak_idx] = max(abs(z));
tau_est = lags(peak_idx) / fs;

%% Estimated Range
R_est = (tau_est * c) / 2;

%% Doppler Estimation via FFT (simplified for multiple chirps)
% Simulate M chirps for Doppler estimation
M = 128;
Y = zeros(M, N);
for m = 1:M
    % Add phase shift due to Doppler
    phase_shift = exp(1j*2*pi*fd*m*T);
    y_m = zeros(1,N);
    y_m(valid_idx) = cos(pi * mu * t_echo(valid_idx).^2) .* real(phase_shift);
    Y(m, :) = y_m;
end

% FFT across chirps at the estimated range bin
range_fft = abs(fft(Y, [], 2));
range_bin = find(max(abs(range_fft'), [], 1) == max(max(abs(range_fft'))), 1);
doppler_fft = fftshift(fft(Y(:, range_bin), M));
[~, doppler_idx] = max(abs(doppler_fft));
fd_est = ((doppler_idx - M/2) / M) * (1/T);  % Doppler freq in Hz

% Estimated Velocity
v_est = fd_est * lambda / 2;

%% Display Results
fprintf("True Range      = %.2f m\n", R_true);
fprintf("Estimated Range = %.2f m\n", R_est);
fprintf("True Velocity   = %.2f m/s\n", v_true);
fprintf("Estimated Velocity = %.2f m/s\n", v_est);

%% Optional: Plot
figure;
subplot(2,1,1);
plot(lags/fs * 1e6, abs(z));
xlabel('Time lag (\mus)'); ylabel('|Correlation|');
title('Matched Filter Output');
grid on;

subplot(2,1,2);
plot((-M/2:M/2-1)*(1/T), abs(doppler_fft));
xlabel('Doppler Frequency (Hz)');
ylabel('|FFT|');
title('Doppler Spectrum');
grid on;

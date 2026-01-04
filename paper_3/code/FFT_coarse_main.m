clear; clc; close all;

% --- Signal parameters ---
fs = 1000;           % sampling frequency (Hz)
N  = 1024;            % number of samples
f0 = 123.45;         % true frequency (Hz)
a  = 1.0;            % amplitude
phi = 0.3;           % phase
w= 2*pi*f0/fs;
vec = w*(0:N-1)';
% --- Simulation parameters ---
%SNR_vec = -10:1:10;   % SNR range in dB
SNR_dB_vec = -10:1:50;
MC_runs    = 10000;     % Monte Carlo trials per SNR

RMSE = zeros(size(SNR_dB_vec));
RMSE_2 = zeros(size(SNR_dB_vec));
CRLB_vals = zeros(size(SNR_dB_vec));

% --- Loop over SNR values ---
for idx = 1:length(SNR_dB_vec)
    SNR_dB = SNR_dB_vec(idx);
    sigma2 = (a^2/2) / (10^(SNR_dB/10));  % noise variance (since cos power = a^2/2)

    f_est_all = zeros(1,MC_runs);
    f_est2_all = zeros(1,MC_runs);
    f_est_all_hat = zeros (1,MC_runs);
    for mc = 1:MC_runs
        n = 0:N-1;
        y = a*cos(2*pi*f0*n/fs + phi) + sqrt(sigma2)*randn(1,N);
        % --- FFT+Parabolic Estimator ---
        [f_est, f_crlb] = freq_est_fft_parabolic(y, fs, a, sigma2);
         f_est2 = freq_est_fft_improved(y, fs);
        f_est_all(mc) = f_est;
        f_est2_all(mc) = f_est2;
    end
    % --- RMSE ---
    RMSE(idx) = sqrt(mean((f_est_all - f0).^2));
    RMSE_2(idx) = sqrt(mean((f_est2_all - f0).^2));
    % --- CRLB (same for all MC runs, so just store one) ---
    CRLB_vals(idx) = sqrt(f_crlb);  % convert to std deviation (Hz)
end
%% --- Plot results ---
figure;
plot(SNR_dB_vec, RMSE, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
plot(SNR_dB_vec, RMSE_2, 'go-', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
plot(SNR_dB_vec, CRLB_vals, 'r--', 'LineWidth', 1.5);
grid on; xlabel('SNR (dB)'); ylabel('Frequency RMSE (Hz)');
legend('FFT+Parabolic RMSE','improve RMSE','CRLB','Location','northeast');
title(sprintf('Frequency Estimation RMSE vs CRLB (N=%d samples)',N));

